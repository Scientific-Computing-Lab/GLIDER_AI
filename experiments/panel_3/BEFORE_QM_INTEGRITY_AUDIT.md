# Panel-III before-QM integrity audit

Audit time: 2026-08-14T19:04:00Z  
Status: **PASS — reference generation authorized**

| Required condition | Evidence | Result |
|---|---|---|
| Preregistration was remotely timestamped before identity selection | commits `5959616459547eff9a9e483043a2aefd953db62a` and documented pre-selection implementation repairs through `b97ca28b5f697cae76708442e89aef8344810330` | PASS |
| Selector consumed only ID and SMILES | identity-only schema `freesolv_id,smiles`; extraction manifest SHA-256 `0bb3359007bd0ff68185f1ad76a9a90ae1fc4cf2f5955041ef176f3ef3a4265d`; identity table SHA-256 `847e45a27e3963cda2ad909b3a381ac7547fd0dd2dbee4cc506ea2f9b5e8abc8` | PASS |
| FreeSolv hydration columns were neither parsed nor copied | extractor reads semicolon fields 0 and 1 only; source manifest records `fields_after_index_1_parsed_or_stored=false` and `hydration_targets_accessed=false` | PASS |
| Exactly 20 identities and 80 fixed geometries | selection manifest: 4 identities in each of 5 strata; geometry manifest: 20 per regime and 80 total | PASS |
| Geometry is the pushed frozen geometry | `configurations.extxyz` SHA-256 `ec62f08d1e8e2ff7ba1b458ef9cf78c5fda957f96437e7c34c854157528152ee`; remote geometry-freeze commit `56bbec13865189a04f81d6960054a1fd21fb443e` | PASS |
| GLIDER checkpoint is exact | SHA-256 `a288afa285128e13cb7a05ba9459f1dce07e79b8638c6e272c2509ae68d7a288` | PASS |
| All GLIDER predictions exist | 80/80 configuration IDs; prediction-tree SHA-256 `1a7c3caa15c99412f5bbcedecef0c7eedad4754337d71e3679e2b73341dd7179` | PASS |
| All eligible comparator predictions exist | MACE-POLAR S/M/L, AIMNet2, MACE-MDP, GFN2-xTB, static-Thole and zero response each cover 80/80; aggregate comparator-tree SHA-256 `751498a255ffc5542a2d26e4e8483251da2b2f8d78e4db8c2fc6ae4d2852f735` | PASS |
| Prediction files use identical fixed probe points | exact array equality against GLIDER for all nine prediction trees | PASS |
| GLIDER numerical audits pass | deterministic rerun maximum discrepancy 0; maximum `|sum Delta q| < 1e-8 e`; rotation audit `pass=true` | PASS |
| No prediction program used a reference target | explicit false fields in all prediction manifests; no reference input accepted by GLIDER inference | PASS |
| No reference artifact exists | no `canonical_panel_3/references` directory and no response-reference or observable-registry file found | PASS |
| Prediction freeze is remotely timestamped | local and `origin/canonical-panel-3` both `f9a87b1d5a91a45f595ccd1891f59ec8d9445fea` | PASS |
| Frozen all-predictions manifest is immutable | SHA-256 `3dc5641a289f48118e0e372c341286573d4a381a2cf2ea9b7ae2b2b80ec8c73e` | PASS |

The MACE-MDP, AIMNet2 and GFN2 scripts retain a prose chronology sentence from
the earlier comparator-closure campaign that says they were executed after
labels existed. Their explicit target-access fields are false. For Panel III,
the Git chronology and this audit establish that execution occurred before any
reference directory or response label existed; the inherited sentence is not
used as Panel-III chronology evidence.

During this audit the manifest-generation command was inadvertently invoked a
second time. It changed only creation-time and parent-commit text in the three
manifest files; it did not run a model or touch a prediction. Those uncommitted
manifest changes were immediately restored byte-for-byte from remote commit
`f9a87b1` before this audit was written. The authoritative all-predictions
manifest hash remains `3dc5641a...8c73e`.

No response-QM calculation may begin until this audit itself is committed and
pushed. The subsequent reference manifest must prove that every reference file
postdates the remote audit commit.
