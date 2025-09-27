# VertexCompositeAnalysis Overview

## Repository Purpose
- Heavy-flavor vertex composite reconstruction and ntupling tailored for Run-3 heavy-ion and pp datasets.
- Provides standalone CMSSW modules orthogonal to HiForest, focused on D-meson and related decay channels (e.g., D⁰→Kπ, D*±→D⁰π, D±→Kππ).

## Package Layout
- `VertexCompositeProducer/`
  - `plugins/`: EDM producers and utilities (e.g., `GenParticleFixer`, `HiHFFilter`, `LumiProducerFromBrilcalc`, `MuonUnpacker`, `NTrackVertexMapper`) plus BuildFile for CMSSW plugin registration.
  - `src/`: Decay-specific candidate builders and fitters (`D0Producer`, `DStarProducer`, `LamC3PProducer`, `V0Producer`, etc.), usually paired with corresponding fitters.
  - `python/`: Producer configuration fragments (general candidate definitions, skims).
  - `test/`: Ready-to-run `cmsRun` configs for 2023 PbPb and 2024 pp campaigns, multi-step tuning (Step0/1/2MVA), CRAB job templates, and database payloads.
- `VertexCompositeAnalyzer/`
  - `plugins/`: Analysis EDM modules (tree/ntuple producers, event-plane utilities, PAT-based analyzers).
  - `python/`: Rich set of analyzer/selectors covering D⁰, D*, D±, Ds, Λc, dimuon, and V0 channels in both tree and ntuple flavors.
  - `test/`: Example analyzer configs per channel (D⁰, JPsi, Λc3P) to drive local validation.
- Top-level `README.md`: setup recipe (e.g., CMSSW_13_2_11) and TODO list for future improvements (MC configs, selection updates, 3-prong optimization).

## Notable Resources
- CRAB configs under `VertexCompositeProducer/test` for data/MC submissions (PbPb2023, pp2024) with multi-crab helper scripts.
- Database payload snapshots (`HeavyIonRPRcd_*.db`, `CentralityTable_*.db`) necessary for heavy-ion workflows.
- Macros and TMVA training materials in `VertexCompositeAnalyzer/macros/` for downstream analysis.

## Current Build Status (CMSSW_14_1_X branch)
- `scram b -j8` currently fails because the build area cannot write to `tmp/el9_amd64_gcc12/...`; ensure the CMSSW working directory is owned or writable before compiling.
- `VertexCompositeAnalyzer/plugins/PATCompositeTreeProducer3.h` includes missing `ValidationUtility.h`; add the header or adjust the include to complete compilation.

## Usage Notes
- Initialize environment with `cmsenv` inside the chosen CMSSW release before running configs.
- Test locally with small `maxEvents` before launching large-scale or CRAB jobs.
