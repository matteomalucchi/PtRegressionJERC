"""PocketCoffea processor for the column-based MC truth derivation.

This is the slim replacement of :mod:`mc_truth_ptreg_jerc.workflow`. Compared to
the old processor:

* it reads **no environment variables** — everything is steered through
  ``workflow_options`` in the configuration file;
* it always computes **both** the ParticleNet and the UParT regression, with
  and without neutrinos, in the same job;
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
    DeltaR_matching : float
        Cone used to match reco jets to gen jets.
    GenJetPtCut : float
        Minimum gen jet pT (applied to the neutrino-inclusive gen jets as well).
    n_jets : int
        Number of leading (raw-pT ordered) reco jets kept per event.
    com_energy : float
        Centre-of-mass energy in TeV, used to reject unphysical jets.
    pnet, upart : bool
        Compute the ParticleNet / UParT regressed response.
    neutrinos : bool
        Also build the neutrino-inclusive gen jet collection and its responses.
    mc_truth_corr : dict or None
        L2Relative correction re-applied on top of the raw jet pT to define the
        ``JEC`` response. When ``None`` the jet pT calibrated by PocketCoffea is
        used instead.
    closure_corrections : dict
        Mapping ``{"PNetReg": corr, "PNetRegNeutrino": corr, "UparTReg": corr,
        "UparTRegNeutrino": corr}`` of freshly derived MC truth corrections to
        apply on top of each regression (closure test). Missing or ``None``
        entries are simply not applied.
    regression_value_when_invalid : {"original", "raw"} or float
        Value used when a regression factor is not positive.
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
        self.DeltaR_matching = opts.get("DeltaR_matching", 0.2)
        self.GenJetPtCut = opts.get("GenJetPtCut", 0)
        self.n_jets = opts.get("n_jets", 3)
        self.com_energy = opts.get("com_energy", 13.6)
        self.pnet = opts.get("pnet", True)
        self.upart = opts.get("upart", True)
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

    def _replace_invalid_regression(
        self, collection_name, pt_field, response_field, invalid_mask
    ):
        """Replace the regressed pT/response wherever *invalid_mask* is True."""
        if self.regression_value_when_invalid == "original":
            return

        collection = self.events[collection_name]
        if self.regression_value_when_invalid == "raw":
            pt_repl = collection["JetPtRaw"]
            response_repl = collection["ResponseRaw"]
        else:
            value = float(self.regression_value_when_invalid)
            pt_repl = ak.full_like(collection[pt_field], value)
            response_repl = ak.full_like(collection[response_field], value)

        self.events[collection_name] = ak.with_field(
            self.events[collection_name],
            ak.where(invalid_mask, response_repl, collection[response_field]),
            response_field,
        )
        self.events[collection_name] = ak.with_field(
            self.events[collection_name],
            ak.where(invalid_mask, pt_repl, collection[pt_field]),
            pt_field,
        )

    def _add_regression(
        self, collection_name, gen_pt, corr_factor, suffix, invalid_mask=None
    ):
        """Add ``JetPt<suffix>`` / ``Response<suffix>`` to *collection_name*.

        *corr_factor* is the multiplicative regression factor applied to the raw
        jet pT. Invalid factors (``invalid_mask``, by default a non positive
        *corr_factor*) and jets outside the training region are replaced
        according to the workflow options.
        """
        pt_field = f"JetPt{suffix}"
        response_field = f"Response{suffix}"

        collection = self.events[collection_name]
        self.events[collection_name] = ak.with_field(
            collection, collection.JetPtRaw * corr_factor, pt_field
        )
        self.events[collection_name] = ak.with_field(
            self.events[collection_name],
            self.events[collection_name][pt_field] / gen_pt,
            response_field,
        )

        if invalid_mask is None:
            invalid_mask = corr_factor <= 0
        self._replace_invalid_regression(
            collection_name, pt_field, response_field, invalid_mask
        )

        if self.invalid_regression_eta is not None:
            outside = (
                abs(self.events[collection_name].RecoEta)
                > self.invalid_regression_eta
            )
            self.events[collection_name] = ak.with_field(
                self.events[collection_name],
                ak.where(
                    outside,
                    self.events[collection_name].ResponseRaw,
                    self.events[collection_name][response_field],
                ),
                response_field,
            )
            self.events[collection_name] = ak.with_field(
                self.events[collection_name],
                ak.where(
                    outside,
                    self.events[collection_name].JetPtRaw,
                    self.events[collection_name][pt_field],
                ),
                pt_field,
            )

    def _apply_closure_correction(self, collection_name, suffix, gen_pt, upart=False):
        """Apply a freshly derived MC truth correction on top of a regression."""
        corr_info = self.closure_corrections.get(suffix)
        if not corr_info:
            return

        pt_field = f"JetPt{suffix}"
        response_field = f"Response{suffix}"
        collection = self.events[collection_name]

        corr = get_mc_truth_corr(
            corr_info,
            collection.RecoEta,
            collection.RecoPhi,
            collection[pt_field],
            pnetreg=not upart,
            upartreg=upart,
        )
        self.events[collection_name] = ak.with_field(
            collection, collection[pt_field] * corr, pt_field
        )
        self.events[collection_name] = ak.with_field(
            self.events[collection_name],
            self.events[collection_name][pt_field] / gen_pt,
            response_field,
        )

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
            # use the jet pT calibrated by PocketCoffea
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
    def apply_object_preselection(self, variation):
        if not self._isMC:
            return

        if "pt_raw" not in self.events["Jet"].fields:
            self.events["Jet"] = ak.with_field(
                self.events["Jet"],
                self.events["Jet"].pt * (1 - self.events["Jet"].rawFactor),
                "pt_raw",
            )

        self.events["JetGood"] = ak.with_field(
            self.events.Jet, self.events.Jet.pt_raw, "ptRaw"
        )
        self.events["JetGood"] = self.events.JetGood[
            ak.argsort(self.events.JetGood.ptRaw, ascending=False)
        ][:, : self.n_jets]

        # reject jets whose energy exceeds the beam energy
        physical_jet_mask = self.events.JetGood.ptRaw * np.cosh(
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

        # --- regular gen jets -----------------------------------------
        gen_matched, reco_matched, _ = object_matching(
            self.events["GenJetGood"], self.events["JetGood"], self.DeltaR_matching
        )
        gen_matched = gen_matched[~ak.is_none(gen_matched, axis=1)]
        reco_matched = reco_matched[~ak.is_none(reco_matched, axis=1)]

        self._build_matched_collection("MatchedJets", gen_matched, reco_matched)

        if self.pnet:
            self._add_regression(
                "MatchedJets",
                gen_matched.pt,
                reco_matched.PNetRegPtRawCorr,
                "PNetReg",
            )
            self._apply_closure_correction("MatchedJets", "PNetReg", gen_matched.pt)

        if self.upart:
            self._add_regression(
                "MatchedJets",
                gen_matched.pt,
                reco_matched.UParTAK4RegPtRawCorr,
                "UparTReg",
            )
            self._apply_closure_correction(
                "MatchedJets", "UparTReg", gen_matched.pt, upart=True
            )

        # --- neutrino-inclusive gen jets -------------------------------
        if not self.neutrinos:
            return

        gen_nu_matched, reco_nu_matched, _ = object_matching(
            self.events["GenJetNeutrino"], self.events["JetGood"], self.DeltaR_matching
        )
        gen_nu_matched = gen_nu_matched[~ak.is_none(gen_nu_matched, axis=1)]
        reco_nu_matched = reco_nu_matched[~ak.is_none(reco_nu_matched, axis=1)]

        self._build_matched_collection(
            "MatchedJetsNeutrino", gen_nu_matched, reco_nu_matched
        )

        if self.pnet:
            # the neutrino regression factor is applied on top of the standard one
            self._add_regression(
                "MatchedJetsNeutrino",
                gen_nu_matched.pt,
                reco_nu_matched.PNetRegPtRawCorr
                * reco_nu_matched.PNetRegPtRawCorrNeutrino,
                "PNetRegNeutrino",
                invalid_mask=(reco_nu_matched.PNetRegPtRawCorr <= 0)
                | (reco_nu_matched.PNetRegPtRawCorrNeutrino <= 0),
            )
            self._apply_closure_correction(
                "MatchedJetsNeutrino", "PNetRegNeutrino", gen_nu_matched.pt
            )

        if self.upart:
            # UParT provides a single neutrino-inclusive factor: unlike PNet the
            # standard regression factor must NOT be multiplied in here.
            self._add_regression(
                "MatchedJetsNeutrino",
                gen_nu_matched.pt,
                reco_nu_matched.UParTAK4RegPtRawCorrNeutrino,
                "UparTRegNeutrino",
                invalid_mask=(reco_nu_matched.UParTAK4RegPtRawCorr <= 0)
                | (reco_nu_matched.UParTAK4RegPtRawCorrNeutrino <= 0),
            )
            self._apply_closure_correction(
                "MatchedJetsNeutrino", "UparTRegNeutrino", gen_nu_matched.pt, upart=True
            )

    def count_objects(self, variation):
        pass
