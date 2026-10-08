# Which labels informed which model?

The published checkpoint is the original five-member response head. Its fitting inputs are the 48 configurations in [training](../experiments/training/), comprising 24 three-water and 24 four-water clusters.

1. The original development study used leave-one-chemistry-out folds to select the response-head construction and training schedule.
2. That checkpoint was frozen before its Panel-I predictions and response references. Seven Panel-I identities had been designed in an earlier sealed campaign. Five further identities were chosen using chemistry-space diversity, without response errors.
3. After Panel I was scored, its 48 labels were combined with the original 48 for a **different later global-prior model family**. That later model failed to improve on Panel II. Its predictions remain under `panel_2/predictions/later_global_prior`.
4. The original GLIDER checkpoint was unchanged for Panels II and III, the non-water geometry challenge and the liquid/size/coupling extensions. Historical identities and results informed the broader research sequence, so the release does not claim that all later analysis decisions were blind to all earlier labels. The non-water contact audit was performed after scoring and shows that most fixed replacements were repulsive; it does not alter the original panel or establish transfer at attractive contacts.
5. The new global-branch and fragment-separation tests were run during author review on 4 October 2026. They are explicitly post hoc and do not replace any original result.
6. After the non-water contact audit, a separate [36-case contact follow-up](../experiments/nonwater_contact/) retained all 12 original solutes and three neighbour species. Geometry was selected without new QM references or response scores. The geometry file and both complete prediction sets were hashed in the [freeze record](../experiments/nonwater_contact/prediction_freeze.json) at 22:07:44 UTC on 7 October 2026, before new response QM was acquired. The checkpoint was unchanged. This records the sequence used for the new scores, but the follow-up is post-critique and is not an independently preregistered replication.

The historical string `OLD_PROSPECTIVE_NOW_DEVELOPMENT` described Panel I's subsequent reuse. It does not mean those labels trained the published GLIDER head. The [later pool](../provenance/later_development_pool/) and [failed later-model records](../provenance/negative_results/) are separated from original training.

## How the identities and geometries were selected

Panel I combines seven designed neutral polyfunctional solutes from a pre-existing sealed list with five max–min diversity selections. The latter use a 0.6 Morgan / 0.4 scaled-Mordred distance with development identities and the seven sealed identities as anchors. Four strata plus one global selection supply the five new identities. The [source selection code and records](../provenance/sampling/) preserve these details.

Panels II and III use their own recorded selection procedures. Panel III accesses FreeSolv identities and connectivity without hydration values. [Original Panel-III selection records](../provenance/panel_3/).

The six liquid-bridge solutes are absent from the original 14-solute response-training set. Four use fully coupled endpoints from earlier replica-exchange simulations and two use newly sampled trajectories. These identities were called development solutes in a separate energy-model campaign. We therefore use the narrower, verifiable statement **outside the original response-training set**, rather than claiming that no earlier research decision ever involved them.

## Historical paths

Paths beginning `/dev/shm/` or `/home/galoren/` in preserved manifests describe where a past calculation ran. They are not instructions for finding files in this release. Active entry points use `experiments/` and repository-relative paths. [Release path map](../provenance/release_path_map.json) records the migration, and [RELEASE_MANIFEST.json](../RELEASE_MANIFEST.json) checks current files. We retain the original manifests unchanged so that the chronology remains auditable.
