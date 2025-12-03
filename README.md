# Hackathon2025
Repository for all contents concerning the 2025 Hackathon by Bayern Innovativ

## Synthetic Graph
The biological background part of the graph (neo4j version 5.25.1) can be downloaded from [LRZ Sync&Share](https://syncandshare.lrz.de/getlink/fiUyxWQz72XE4ttzusb9qn/biological_background.dump).

The `SYNTETIC_singleScript.cypher`file produces synthetic patient data with the same connection and attribute types as the real patient data.

### Online versions
Online versions of the synthetic graph can be found in two places:
- [Kamatera Server](http://83.229.84.12:7474/browser/) User: tumaiReadonly PW: MAKEATHON2024 (_for exploration only_)
- [Azure Server](http://52.157.242.238:7474/browser/) User: neo4j PW: Neo4JJJJJ (_get access to the azure vm for code execution via the slack_)

## About the Graph
It is a highly modified version of the [clinical knowledge graph](https://ckg.readthedocs.io/en/latest/) and has not been published yet.
As it contains proprietary ontologies such as SNOMED-CT, the graph is explicitly __NOT OPEN SOURCE__.

The documentation for the graph can be found in `graph_README.md`. 

## Featurecloud
We use the [FeatureCloud](featurecloud.ai) technology to grant federated access to the patient data. A tutorial on how to make your script compatible with FeatureCloud can be found [here](https://github.com/careforrare/PersonalizedMedicine)

## Submission
Final submission of the Hackathon results will happen via [FTAPI](https://datatransfer.bayern-innovativ.de/submit/lara_kronester) a guid on how to do this can be found in the repo.
__All results have to be handed in 12:00 on December 6th__

## Contact
Join the [Slack](https://join.slack.com/t/hackathon2025gruppe/shared_invite/zt-3k6cgyvgy-ccsf6fbIT6Yth4g2UqHtkw).
Should you have any questions, contact [Henrik Otterstedt](mailto:Henrik.Otterstedt@med.uni-muenchen.de) or open an issue.
