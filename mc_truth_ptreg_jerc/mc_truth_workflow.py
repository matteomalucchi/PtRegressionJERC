"""PocketCoffea processor for the column-based MC truth derivation.

This is the slim replacement of :mod:`mc_truth_ptreg_jerc.workflow`. Compared to
the old processor:

* it reads **no environment variables** -- everything is steered through
  ``workflow_options`` in the configuration file;
* the pT regressions are **not** computed here: they are applied by the
  PocketCoffea ``JetsCalibrator``, configured by
  ``params/jets_calibration_mc_truth.yaml``, in the same way as in the
  ``met_ptreg_performance`` configuration. This processor only clones the
  ``Jet`` collection once per regression before the calibrators run, and then
  reads back the regressed pT;
* both ParticleNet and UParT are processed in the same job, with and without
  neutrinos;
* it does **not** create per-eta / per-pT collections or categories: the
  binning is done offline on the dumped columns, so a single ``baseline``
  category is enough;
* the parton/hadron flavour of the matched gen jet is kept as a column so that
  the flavour splitting can be done at plotting time.

The output is dumped as parquet files (one per chunk) via the
``dump_columns_as_arrays_per_chunk`` workflow option, exactly like the JER MC
configuration does.
"""

import awkward as ak
import numpy as np

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.deltaR_matching import object_matching, deltaR_matching_nonunique

from mc_truth_ptreg_jerc.custom_functions import get_mc_truth_corr


class MCTruthProcessor(BaseProcessorABC):
    """Match reco jets to gen jets and compute the response of every correction.

    Workflow options
    ----------------
    regressions : list of dict
        One entry per pT regression, e.g.::

            {"name": "PNetReg", "collection": "JetPNetReg", "neutrinos": False}

        ``collection`` must be the collection that
        ``params/jets_calibration_mc_truth.yaml`` maps to the corresponding
        regression jet type: the processor clones ``Jet`` into it before the
        calibrators run, and the ``JetsCalibrator`` replaces its pT with the
        regressed one. ``name`` is the suffix of the output columns
        (``JetPt<name>`` / ``Response<name>``) and ``neutrinos`` says whether the
        response is measured against the neutrino-inclusive gen jets.
    DeltaR_matching : float
        Cone used to match reco jets to gen jets.
    GenJetPtCut : float
        Minimum gen jet pT (applied to the neutrino-inclusive gen jets as well).
    n_jets : int
        Number of leading (raw-pT ordered) reco jets kept per event.
    com_energy : float
        Centre-of-mass energy in TeV, used to reject unphysical jets.
    neutrinos : bool
        Also build the neutrino-inclusive gen jet collection and its responses.
    mc_truth_corr : dict or None
        L2Relative correction re-applied on top of the raw jet pT to define the
        ``JEC`` response. When ``None`` the jet pT calibrated by PocketCoffea is
        used instead.
    closure_corrections : dict
        Mapping ``{"<regression name>": corr}`` of freshly derived MC truth
        corrections to apply on top of each regression (closure test). Missing
        or ``None`` entries are simply not applied.
    regression_value_when_invalid : {"original", "raw"} or float
        Value used when the regressed pT is not positive, which is how the
        calibrator leaves the jets whose regression factor is not valid.
    invalid_regression_eta : float or None
        Above this ``|eta|`` the regression output is replaced by the raw jet
        (the taggers are not trained there). ``None`` disables the replacement.

    The ``partonFlavour`` / ``hadronFlavour`` fields of the gen jet are carried
    over to the matched collections automatically, so they can simply be listed
    in the ``ColOut`` of the configuration file.
    """

    def __init__(self, cfg: Configurator):
        super().__init__(cfg)

        opts = self.workflow_options
        self.regressions = opts.get("regressions", [])
        self.DeltaR_matching = opts.get("DeltaR_matching", 0.2)
        self.GenJetPtCut = opts.get("GenJetPtCut", 0)
        self.n_jets = opts.get("n_jets", 3)
        self.com_energy = opts.get("com_energy", 13.6)
        self.neutrinos = opts.get("neutrinos", True)
        self.mc_truth_corr = opts.get("mc_truth_corr", None)
        self.closure_corrections = opts.get("closure_corrections", {}) or {}
        self.regression_value_when_invalid = opts.get(
            "regression_value_when_invalid", "raw"
        )
        self.invalid_regression_eta = opts.get("invalid_regression_eta", 4.7)

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _regressions_for(self, neutrinos):
        """The regressions measured against the (neutrino-inclusive) gen jets."""
        return [
            regression
            for regression in self.regressions
            if bool(regression.get("neutrinos", False)) == neutrinos
        ]

    def add_neutrinos_to_genjets(self, genjets, neutrinos):
        """Add the 4-momentum of the matched gen neutrinos to the gen jets."""
        neutrinos_matched = deltaR_matching_nonunique(genjets, neutrinos, 0.4)

        # sum the 4-vectors of all the neutrinos matched to a given gen jet
        sum_px = ak.sum(neutrinos_matched.px, axis=-1)
        sum_py = ak.sum(neutrinos_matched.py, axis=-1)
        sum_pz = ak.sum(neutrinos_matched.pz, axis=-1)
        sum_e = ak.sum(neutrinos_matched.energy, axis=-1)

        sum_pt = ak.nan_to_num(np.sqrt(sum_px**2 + sum_py**2), nan=0)
        sum_eta = ak.nan_to_num(
            np.arctanh(sum_pz / np.sqrt(sum_pt**2 + sum_pz**2)), nan=0
        )
        sum_phi = ak.nan_to_num(np.arctan2(sum_py, sum_px), nan=0)
        sum_mass = ak.nan_to_num(
            np.sqrt(sum_e**2 - sum_pt**2 - sum_pz**2), nan=0
        )

        neutrinos_sum = ak.zip(
            {"pt": sum_pt, "eta": sum_eta, "phi": sum_phi, "mass": sum_mass},
            with_name="PtEtaPhiMLorentzVector",
        )

        genjets_with_neutrinos = genjets + neutrinos_sum

        # keep all the original gen jet fields (flavour included) and only
        # overwrite the kinematics
        out = genjets
        for field in ("pt", "eta", "phi", "mass"):
            out = ak.with_field(out, getattr(genjets_with_neutrinos, field), field)
        return out

    def _regressed_pt_field(self, name):
        """Name of the field holding the regressed pT on the Jet collection."""
        return f"pt{name}"

    def _add_regressed_response(self, collection_name, name, regressed_pt, gen_pt):
        """Add ``JetPt<name>`` / ``Response<name>`` from an already regressed pT.

        The regression itself was applied by the ``JetsCalibrator``; the only
        thing left to do here is to recover the analysis-level fallbacks: the
        calibrator leaves ``pt = 0`` where the regression factor is not valid,
        and the taggers are not trained above ``invalid_regression_eta``.

        Note that the validity is judged on the regressed pT rather than on the
        individual regression factors, which the calibrator does not expose.
        The two are equivalent for the invalid value the taggers actually write
        (zero); a pair of *negative* factors multiplied together would not be
        caught, but that does not occur in practice.
        """
        pt_field = f"JetPt{name}"
        response_field = f"Response{name}"
        collection = self.events[collection_name]
        pt = regressed_pt

        if self.regression_value_when_invalid != "original":
            if self.regression_value_when_invalid == "raw":
                fallback = collection.JetPtRaw
            else:
                fallback = ak.full_like(pt, float(self.regression_value_when_invalid))
            pt = ak.where(pt > 0, pt, fallback)

        if self.invalid_regression_eta is not None:
            pt = ak.where(
                abs(collection.RecoEta) > self.invalid_regression_eta,
                collection.JetPtRaw,
                pt,
            )

        collection = ak.with_field(collection, pt, pt_field)
        collection = ak.with_field(collection, pt / gen_pt, response_field)
        self.events[collection_name] = collection

    def _apply_closure_correction(self, collection_name, name, gen_pt):
        """Apply a freshly derived MC truth correction on top of a regression."""
        corr_info = self.closure_corrections.get(name)
        if not corr_info:
            return

        pt_field = f"JetPt{name}"
        response_field = f"Response{name}"
        collection = self.events[collection_name]

        upart = "UparT" in name or "UParT" in name
        corr = get_mc_truth_corr(
            corr_info,
            collection.RecoEta,
            collection.RecoPhi,
            collection[pt_field],
            pnetreg=not upart,
            upartreg=upart,
        )
        collection = ak.with_field(collection, collection[pt_field] * corr, pt_field)
        collection = ak.with_field(
            collection, collection[pt_field] / gen_pt, response_field
        )
        self.events[collection_name] = collection

    def _build_matched_collection(self, collection_name, gen_jets, reco_jets):
        """Create the matched collection with the raw and JEC responses."""
        self.events[collection_name] = ak.with_field(
            gen_jets, reco_jets.eta, "RecoEta"
        )
        self.events[collection_name] = ak.with_field(
            self.events[collection_name], reco_jets.phi, "RecoPhi"
        )
        self.events[collection_name] = ak.with_field(
            self.events[collection_name], reco_jets.pt_raw, "JetPtRaw"
        )
        self.events[collection_name] = ak.with_field(
            self.events[collection_name],
            self.events[collection_name].JetPtRaw / gen_jets.pt,
            "ResponseRaw",
        )

        if self.mc_truth_corr:
            # re-apply the standard L2Relative on top of the raw pT
            corr = get_mc_truth_corr(
                self.mc_truth_corr,
                self.events[collection_name].RecoEta,
                self.events[collection_name].RecoPhi,
                self.events[collection_name].JetPtRaw,
                pnetreg=False,
            )
            jet_pt_jec = self.events[collection_name].JetPtRaw * corr
        else:
            # use the jet pT calibrated by the JetsCalibrator
            jet_pt_jec = reco_jets.pt

        self.events[collection_name] = ak.with_field(
            self.events[collection_name], jet_pt_jec, "JetPtJEC"
        )
        self.events[collection_name] = ak.with_field(
            self.events[collection_name],
            self.events[collection_name].JetPtJEC / gen_jets.pt,
            "ResponseJEC",
        )

    # ------------------------------------------------------------------
    # PocketCoffea interface
    # ------------------------------------------------------------------
    def process_extra_after_skim(self):
        """Prepare the collections that the ``JetsCalibrator`` will regress.

        This runs before the calibrators are initialized. Each regression gets
        its own clone of the raw ``Jet`` collection; the calibrator then
        replaces its pT with the regressed one (see
        ``params/jets_calibration_mc_truth.yaml``).
        """
        if not self._isMC:
            return

        if "pt_raw" not in self.events["Jet"].fields:
            self.events["Jet"] = ak.with_field(
                self.events["Jet"],
                self.events["Jet"].pt * (1 - self.events["Jet"].rawFactor),
                "pt_raw",
            )

        for regression in self.regressions:
            self.events[regression["collection"]] = self.events["Jet"]

    def process_extra_after_calibrators(self, variation):
        """Copy the regressed pT computed by the calibrators onto ``Jet``.

        All the regressed collections are clones of ``Jet`` and none of them is
        re-sorted (``sort_by_pt: False``), so they stay index-aligned and the
        regressed pT can simply travel as extra fields of the jets that are
        matched to the gen jets. This way the matching is done once and all the
        corrections are compared on exactly the same jets.
        """
        if not self._isMC:
            return

        for regression in self.regressions:
            self.events["Jet"] = ak.with_field(
                self.events["Jet"],
                self.events[regression["collection"]].pt,
                self._regressed_pt_field(regression["name"]),
            )

    def apply_object_preselection(self, variation):
        if not self._isMC:
            return

        self.events["JetGood"] = self.events.Jet[
            ak.argsort(self.events.Jet.pt_raw, ascending=False)
        ][:, : self.n_jets]

        # reject jets whose energy exceeds the beam energy
        physical_jet_mask = self.events.JetGood.pt_raw * np.cosh(
            self.events.JetGood.eta
        ) < (self.com_energy * 1000) / 2
        self.events["JetGood"] = self.events.JetGood[physical_jet_mask]

        self.events["GenJetGood"] = self.events.GenJet[
            self.events.GenJet.pt > self.GenJetPtCut
        ]

        if self.neutrinos:
            neutrinos = self.events["GenPart"][
                (abs(self.events.GenPart.pdgId) == 12)
                | (abs(self.events.GenPart.pdgId) == 14)
                | (abs(self.events.GenPart.pdgId) == 16)
            ]
            self.events["GenJetNeutrino"] = self.add_neutrinos_to_genjets(
                self.events.GenJet, neutrinos
            )
            self.events["GenJetNeutrino"] = self.events.GenJetNeutrino[
                self.events.GenJetNeutrino.pt > self.GenJetPtCut
            ]

    def process_extra_after_presel(self, variation) -> ak.Array:
        if not self._isMC:
            return

        for collection_name, gen_collection, neutrinos in (
            ("MatchedJets", "GenJetGood", False),
            ("MatchedJetsNeutrino", "GenJetNeutrino", True),
        ):
            if neutrinos and not self.neutrinos:
                continue

            gen_matched, reco_matched, _ = object_matching(
                self.events[gen_collection],
                self.events["JetGood"],
                self.DeltaR_matching,
            )
            gen_matched = gen_matched[~ak.is_none(gen_matched, axis=1)]
            reco_matched = reco_matched[~ak.is_none(reco_matched, axis=1)]

            self._build_matched_collection(collection_name, gen_matched, reco_matched)

            for regression in self._regressions_for(neutrinos):
                name = regression["name"]
                self._add_regressed_response(
                    collection_name,
                    name,
                    reco_matched[self._regressed_pt_field(name)],
                    gen_matched.pt,
                )
                self._apply_closure_correction(collection_name, name, gen_matched.pt)

    def count_objects(self, variation):
        pass
