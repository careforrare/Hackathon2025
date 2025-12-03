# Explaining the different Node Classes

#### Chromosome

A chromosome is the most high level unit of organisation of DNA. The human body has 46 of them. Each chromosome has a duplicate, meaning that there are 23(24) different types of them. (The 24 different types come from counting the X and Y chromosome separately).

Additionally the graph also includes the mitochondrial chromosome bringing the total to 25.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | Either the number or "X", "Y", "MT" | String |
| id | The same as the name | String |


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Transcript / Gene | LOCATED_IN | Chromosome | start : int <br> end : int <br> strand : String ('+' or '-') |
| Known_variant | VARIANT_FOUND_IN_CHROMOSOME | Chromosome | N/A |


#### Gene

A gene is the basic functional unit of the human genome. It is made up of a sequence of nucleic acids in so called exons and introns. Where the exons are the part of the gene that get transcribed into a protein, whilst the introns are "filler" between the exons.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| family | optional: name of the gene familiy | String
| synonyms | list of alternative gene identifiers e.g ENSEMBL id | List of Strings
| id | Gene name | String
| name | Short description of the gene and if known its function | String


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Gene | LOCATED_IN | Chromosome | start : int <br> end : int <br> strand : String ('+' or '-') |
| Gene | ASSOCIATED_WITH | Disease | number_publications : String |
| Gene | TRANSCRIBED_INTO | Transcript | N/A |
| Gene | TRANSLATED_INTO | Protein | N/A |
| Known_variant | VARIANT_FOUND_IN_GENE | Gene | N/A
| Biological_sample | HAS_DAMAGE | Gene | cadd : String
| Phenotype | ASSOCIATED_WITH | Gene | disease: String (OMIM or ORPHA code of diasease that supports the association |


#### Transcript

The transcript is the part of the gene that gets turned into a protein. In a process called splicing the introns get removed from the genetic sequence and we are only left with the exon sequences to make up the transcript. For a given gene there can be different transcripts depending on which order the exons get added together/if they get added at all. (Alternative Splicing)

A transcript then gets translated into a protein. Here the nucleic acid triplets each code for one of twenty different amino acid that make up proteins.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| assembly | Which assembly of the human genome the transcript stems from | String
| name | Short description of the transcript i.e. what sort of molecule does it code for | String
| id | the transcript Genbank id | String
| class | which class of transcript is it (i.e. what type of RNA) e.g. lncRNA, protein_coding, with_protein | String


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Gene | TRANSCRIBED_INTO | Transcript | N/A |
| Transcript | LOCATED_IN | Chromosome | start : int <br> end : int <br> strand : String ('+' or '-') |
| Transcript | TRANSLATED_INTO | Protein | N/A | 


#### Protein

Proteins are the molecules that fulfill nearly all functionalities needed in the body. They are made up of 20 different amino acids. A proteins function is highly dependent on its structure, which it folds into due to the physical interactions between the different amino acids.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | Protein name | String
| id | UniProt id | String
| sequence | Amino Acid sequence of the protein | String
| synonyms | list of alternative names for the protein e.g. PDB, Ensembl | List of Strings
| accession | Gene name that it is translated from | String


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Protein | DETECTED_IN_PATHOLOGY_SAMPLE | Disease | not_detected : String <br> linkout : String (Hyperlink to the source of the detection) <br> negative_prognosis_logrank_pvalue : float <br> source : String <br> positive_prognosis_logrank_pvalue : float <br> expression_low : float <br> expression_medium : float <br> expression_high : float |
| Protein | ASSOCIATED_WITH | Cellular_component | evidence_type : String <br> source : String <br> score : float |
| Protein | ASSOCIATED_WITH | Tissue | evidence_type : String <br> source : String <br> score : float |
| Protein | ASSOCIATED_WITH | Disease | evidence_type : String <br> source : String <br> score : float |
| Protein | ASSOCIATED_WITH | Biological_process | evidence_type : String <br> source : String <br> score : float |
| Protein | ASSOCIATED_WITH | Molecular_function | evidence_type : String <br> source : String <br> score : float |
| Protein | COMPILED_INTERACTS_WITH | Protein | score : float <br> source : String <br> evidence : String <br> scores : List of String |
| Protein | ACTS_ON | Protein | directionallity : String (possible values: true/false) <br> action : String <br> source : String <br> score : float |
| Protein | CURATED_INTERACTS_WITH | Protein | interaction_type : String <br> source : String <br> method : String <br> evidence : List of Strings <br> score : float |
| Protein | MENTIONED_IN_PUBLICATION | Publication | N/A |
| Protein | HAS_MODIFIED_SITE | Modified_protein | N/A |
| Protein | ANNOTATED_IN_PATHWAY | Pathway |  cellular_component : String <br> organism : String <br> source : String <br> evidence : String |
| Protein | IS_SUBUNIT_OF | Complex | source : String <br> cell_lines : String <br> publication : String <br> evidences : List of Strings |
| Protein | IS_BIOMARKER_OF_DISEASE | Disease | assay : String <br> sex : String <br> units : String <br> age_range : String (of type 'start_age-end_age') <br> age_units : String <br> notes : String <br> source : String <br> is_routine : String (possible values: true/false) <br> normal_range : String (of type 'start-end') <br> is_used_in_clinic : String (possible values: true/false) <br> reference : String |
| Protein | IS_QCMARKER_IN_TISSUE | Tissue | class : String |
| Peptide | BELONGS_TO_PROTEIN | Protein | source : String |
| Transcript | TRANSLATED_INTO | Protein | N/A | 
| Gene | TRANSLATED_INTO | Protein | N/A |
| Known_variant | VARIANT_FOUND_IN_PROTEIN | Protein | N/A |
| Known_variant | CURATED_AFFECTS_INTERACTION_WITH | Protein | effect : String <br> evidence : String <br> interaction : String <br> internal_id : String <br> source : String |
| Biological_sample | HAS_QUANTIFIED_PROTEIN | Protein | source : String (possible values: urine/blood) <br> score : float (log10 quantification) |
| Functional_region | FOUND_IN_Protein | Protein | start : int (local position in protein sequence) <br> end : int (local position in protein sequence) <br> source : String <br> alignment : String |
| Modified_protein | IS_SUBSTRATE_OF | Protein | evidence_type : String <br> regulation : String <br> score : float <br> source : String |
| Metabolite | ASSOCIATED_WITH | Protein | N/A |


#### Functional_region

The specific part of the protein that makes it accomplish its function.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | Name of the functional region | String |
| description | Short description of the region and what its function is | String |
| id | PFAM accession number | String |
| source | database the region was taken from e.g. PFAM |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Functional_region | FOUND_IN_Protein | Protein | start : int (local position in protein sequence) <br> end : int (local position in protein sequence) <br> source : String <br> alignment : String |

#### Modified_protein

A protein that has a post-translational-modification (PTM) done to itself. These modification are chemical alterations of some of the amino acids and lead to altered protein properties. Examples of modifications include methylation or phosphorilation of certain amino acids.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| id | Id of the protein according to the standart nomenclature for modified proteins <br> e.g Prot_S1234-p | String |
| sequence_window | Part of the protein sequence around the modification | String |
| protein | Name of the protein that was modified | String |
| source | Source that found the modification e.g experimentally_identified | String |
| residue | Modified residue | String |
| position | Position in the protein sequence that was modified | int |


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Protein | HAS_MODIFIED_SITE | Modified_protein | N/A |
| Peptide | HAS_MODIFIED_SITE | Modified_protein | N/A |
| Modified_protein | IS_SUBSTRATE_OF | Protein | evidence_type : String <br> regulation : String <br> score : float <br> source : String |
| Modified_protein | MENTIONED_IN_PUBLICATION | Publication | N/A |
| Modified_protein | HAS_MODIFICATION | Modification | source : String |


#### Modification

Post-translational modification (PTM) of a protein that alters its physio-chemical properties. Hierarchical representation of the different modification that are possible according to the [mod ontology](https://www.ebi.ac.uk/ols4/ontologies/mod).

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | Name of the modification | String |
| synonyms | alternative names for the modification | List of Strings |
| id | mod ontology identifier of the modification | String |
| description | Description of the modification | String |


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Modification | HAS_PARENT | Modification | N/A |
| Modified_protein | HAS_MODIFICATION | Modification | source : String |


#### Peptide

A short sequence of amino acids, usually a fraction of a protein produced via some sort of cleavage event e.g. produced via a metabolic process.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| id | Amino acid sequence of the peptide | String |
| type | type of the peptide i.e. which enzyme did the cleavage <br> e.g. tryptic peptide | String |
| unique | Does this peptide originate only from one protein sequence | String (possible values: True/False) |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Peptide | BELONGS_TO_PROTEIN | Protein | source : String |
| Peptide | HAS_MODIFIED_SITE | Modified_protein | N/A |

#### Complex

A complex of proteins that form a distinct structure to accomplish its function.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | Name of the protein complex | String |
| id | CORUM complex id | String |
| source | database the complex was taken from e.g. CORUM | String |
| organism | unique number identifier for the organism in which the complex exists <br> eg. 9606 for human | String | 
| synonyms | List of alternative names for the complex | List of Strings |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Protein | IS_SUBUNIT_OF | Complex | source : String <br> cell_lines : String <br> publication : String <br> evidences : List of Strings |
| Complex | ASSOCIATED_WITH | Biological_process | evidence_type : String <br> source : String <br> score : float

#### Known_variant

A variant of a protein. The variation is due to a mutation of the genes nucleic acid sequence. The different variants can be classified by exploring the effect that they have on the protein. These different effects are encoded in the `effect` attribute of the known_variant nodes.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| pvariant_id | protein_id.amino_acid_change | String |
| external_id | id the variant has in an external database such as TOPmed | String |
| original_source | List of databases the variant is annotated in besides UniProt | List of Strings |
| source | database the variant is taken from | String |
| effect | Effect the variant has e.g. missense variant | String |
| clinical_relevance | clinical_relevance of the mutation <br> if there is none "-" | String |
| disease | Disease the variant is associated with <br> if there is none "-" | String |
| id | notation of the mutation i.e. genomic_location+wt_nucleic_acid'>'mutated_nucleic_acid | String |
| alternative_names | List of alternative names for the Variant | List of String |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Known_variant | VARIANT_FOUND_IN_GENE | Gene | N/A |
| Known_variant | VARIANT_FOUND_IN_PROTEIN | Protein | N/A |
| Known_variant | CURATED_AFFECTS_INTERACTION_WITH | Protein | effect : String <br> evidence : String <br> interaction : String <br> internal_id : String <br> source : String |
| Known_variant | VARIANT_FOUND_IN_CHROMOSOME | Chromosome | N/A |
| Known_variant | VARIANT_FOUND_IN_GWAS | GWAS_study | pvalue : String (in scientific notation i.e. 1E-26) <br> frequency : float (if not set String "NR") <br> source : String <br> trait : String <br> odds_ratio : String (if not set "NR" otherwise the string of a float) |
| Known_variant | VARIANT_IS_CLINICALLY_RELEVANT | Clinically_relevant_variant | source : String |


#### Clinically_relevant_variant

A mutation that is clinically relevant in so far, as it causes a change in the transcribed protein which is associated with a disease, either as a marker or as the (suspected) cause.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| reference | wildtype nucleic acid | String |
| chromosome | chromosome the variant is on | String |
| alternative | new amino acid | String |
| position | position on the chromosome the changed nucleic acid is at | int |
| source | database the variant is taken from | String |
| id | protein_id.amino_acid_change | String |
| alternative_names | List of alternative names for the variant | List of String |
| oncogeneicity | does this variant cause cancer | String |
| effect | effect the variant has | String |


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Known_variant | VARIANT_IS_CLINICALLY_RELEVANT | Clinically_relevant_variant | source : String |
| Clinically_relevant_variant | ASSOCIATED_WITH | Disease | evidence_type : String <br> source : String <br> evidence : String (curated vs not) <br> number_publications : int |


#### Disease

A disease that a patient might be diagnosed with. Ranging from high-level diagnoses to very specific diseases. e.g `shrimp allergy` is less specific than `tiger prawn allergy`. Taken from the [disease ontology](https://www.disease-ontology.org/).

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | name of the disease | String |
| description | short description of the disease | String |
| id | DOID of the disease | String |
| synonyms | List of alternative names and identifiers for the disease | List of String |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Disease | HAS_PARENT | Disease | N/A |
| Gene | ASSOCIATED_WITH | Disease | number_publications : String |
| Protein | ASSOCIATED_WITH | Disease | evidence_type : String <br> source : String <br> score : float |
| Protein | DETECTED_IN_PATHOLOGY_SAMPLE | Disease | not_detected : String <br> linkout : String (Hyperlink to the source of the detection) <br> negative_prognosis_logrank_pvalue : float <br> source : String <br> positive_prognosis_logrank_pvalue : float <br> expression_low : float <br> expression_medium : float <br> expression_high : float |
| Protein | IS_BIOMARKER_OF_DISEASE | Disease | assay : String <br> sex : String <br> units : String <br> age_range : String (of type 'start_age-end_age') <br> age_units : String <br> notes : String <br> source : String <br> is_routine : String (possible values: true/false) <br> normal_range : String (of type 'start-end') <br> is_used_in_clinic : String (possible values: true/false) <br> reference : String |
| Clinically_relevant_variant | ASSOCIATED_WITH | Disease | evidence_type : String <br> source : String <br> score : String (curated vs not) <br> number_publications : String |
| Phenotype | MAPS_TO | Disease | N/A |
| Biological_sample | HAS_DISEASE | Disease | genetically_confirmed : String (possible values: True/False) |
| Disease | MAPS_TO | Clinical_variable | N/A |
| Disease | MENTIONED_IN_PUBLICATION | Publication | N/A |
| Metabolite | ASSOCIATED_WITH | Disease | N/A |

#### Phenotype

In general a phenotype is the outwardly visible expression of the "genetic configuration" of an individual. The phenotypes in the graph are taken from the *[Human Phenotype Ontology](https://hpo.jax.org/).* They provide a standardized hierarchical vocabulary of phenotypic abnormalities encountered in the context of human diseases.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | name of the phenotype | String |
| synonyms | list of alternative names and identifiers for the phenotype | List of String |
| id | HPO id | String |
| description | Short Description of the phenotype | String |


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Phenotype | HAS_PARENT | Phenotype | N/A |
| Phenotype | MAPS_TO | Disease | N/A |
| Phenotype | MAPS_TO | Clinical_variable | N/A |
| Biological_sample | HAS_PHENOTYPE | Phenotype | N/A |
| Phenotype | ASSOCIATED_WITH | Gene | disease: String (OMIM or ORPHA code of diasease that supports the association |


#### Tissue

Description of the tissue (i.e. an assembly of similar cells and their extracellular matrix) according to the [BRENDA tissue ontology](https://en.wikipedia.org/wiki/BRENDA_tissue_ontology).

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| id | BTO id | String |
| description | short description of the tissue | String |
| name | name of the tissue | String |
| synonyms | list of alternative names and identifiers | List of Strings |


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Tissue | HAS_PARENT | Tissue | N/A |
| Protein | ASSOCIATED_WITH | Tissue | evidence_type : String <br> source : String <br> score : float |
| Protein | IS_QCMARKER_IN_TISSUE | Tissue | class : String |
| Tissue | MENTIONED_IN_PUBLICATION | Publication | N/A |


#### Biological_process

The `Biological Process` part of the [gene ontology](geneontology.org). It is a hierarchical classification of larger biological processes that are accomplished by multiple molecular activity. It is important to note that a biological process is not equivalent to a pathway. This is due to the fact that two different pathways can both arrive at the same outcome via different mechanisms.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | Name of the process | String |
| description | Short Description of the biological process | String |
| id | Gene Ontology id | String |
| synonyms | List of alternative names and identifiers | List of Strings |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Biological_process | HAS_PARENT | Biological_process | N/A |
| Protein | ASSOCIATED_WITH | Biological_process | evidence_type : String <br> source : String <br> score : float |
| Complex | ASSOCIATED_WITH | Biological_process | evidence_type : String <br> source : String <br> score : float |

#### Cellular_component

The `Cellular Component` part of the [gene ontology](geneontology.org). It is a hierarchical classification of the location relative to cellular compartments, which a macromolecular machine 'inhabits'.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | name of the GO Term | String |
| description | short description | String |
| id | Gene Ontology id | String |
| synonyms | List of alternative names and identifiers | List of Strings |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Protein | ASSOCIATED_WITH | Cellular_component | evidence_type : String <br> source : String <br> score : float |
| Cellular_component | HAS_PARENT | Cellular_component | N/A | 
| Cellular_component | MENTIONED_IN_PUBLICATION | N/A |

#### Molecular_function

The `Molecular Function` part of the [gene ontology](geneontology.org). It is a hierarchical classification of the molecular-level activities a gene product performs.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | name of the GO Term | String |
| description | short description | String |
| id | Gene Ontology id | String |
| synonyms | List of alternative names and identifiers | List of Strings |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Molecular_function | HAS_PARENT | Molecular_function | N/A |
| Protein | ASSOCIATED_WITH | Molecular_function | evidence_type : String <br> source : String <br> score : float |

#### GWAS_study

A genome-wide association study (GWAS) that observes the association between genetic variants and certain traits, such as diseases.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| date | Date the study was published in "year-month-day" format | String |
| sample_size | short description of the population in the study | String |
| replication_size | size and description of the group used for validation | String |
| trait | Trait studied by the GWAS | String |
| id | unique id  of the study | String |
| title | Title of the study | String |


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Known_variant | VARIANT_FOUND_IN_GWAS | GWAS_study | pvalue : String (in scientific notation i.e. 1E-26) <br> frequency : float (if not set String "NR") <br> source : String <br> trait : String <br> odds_ratio : String (if not set "NR" otherwise the string of a float) |
| GWAS_study | PUBLISHED_IN | Publication | N/A |

#### Metabolite

A small molecule that is an intermediate or end product of metabolism.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| id | Human Metabolite Database id | String |
| name | Name of the metabolite | String |
| sub_class | specific sub class of molecules the metabolite belongs to | String |
| synonyms | List of alternative names and identifiers | List of String |
| direct_parent | precursor  molecule the metabolite is derived from | String |
| super_class | super class of metabolites it belongs to | String |
| description | short description of the molecule | String |
| chemical_fromula | chemical formula of the metabolite | String |
| average_molecular_weight | average molecular weight of the metabolite | Float | 
| class | class of metabolites it belongs to | String |
| kingdom | kingdom of metabolites it belongs to | String | 


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Metabolite | ANNOTATED_IN_PATHWAY | Pathway | cellular_component : String <br> organism : String <br> evidence : String <br> source : String | 
| Metabolite | ASSOCIATED_WITH | Disease | N/A |
| Metabolite | ASSOCIATED_WITH | Protein | N/A |


#### Pathway

A series of actions of different molecules in a cell that leads to a certain product or a change in the cell. This can be the assembly of new molecules, such as fats, turning specific genes on and off, or encourage the cell to move.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | name of the pathway | String |
| description | description of the pathways function | String |
| linkout | hyperlink to the reactome entry for the pathway | String |
| id | reactome id of the pathway | String |
| source | database the pathway was taken from | String |
| organism | organism the pathway was detected in | String |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Metabolite | ANNOTATED_IN_PATHWAY | Pathway | cellular_component : String <br> organism : String <br> evidence : String <br> source : String | 
| Protein | ANNOTATED_IN_PATHWAY | Pathway |  cellular_component : String <br> organism : String <br> source : String <br> evidence : String |

#### Publication

A node representing a distinct publication.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| id | unique number assigned to each publication | String |
| journal | name of the publications journal | String |
| volume | volume of the journal the publication is in | String |
| issue | issue of the journal the publication is in | String |
| DOI | Digital Object Identifier of the publication | String |
| page | Page in the journal the publication starts | String |
| linkout | hyperlink to the publication | String |
| year | year the publication was published in | String |
| PMC_id | PMC id of the publication | String |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Cellular_component | MENTIONED_IN_PUBLICATION | Publication | N/A |
| Protein | MENTIONED_IN_PUBLICATION | Publication | N/A |
| Tissue | MENTIONED_IN_PUBLICATION | Publication | N/A |
| Disease | MENTIONED_IN_PUBLICATION | Publication | N/A |
| Functional_region | MENTIONED_IN_PUBLICATION | Publication | soruce : String |
| Modified_protein | MENTIONED_IN_PUBLICATION | Publication | N/A |
| GWAS_study | PUBLISHED_IN | Publication | N/A |


#### Clinical_variable

The clinical variable node class is the graph representation of [SNOMED-CT](https://www.snomed.org/), which is a licensed ontology of standardized medical language describing both procedures and diagnoses.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| name | Name of the SNOMED term | String |
| id | SNOMED id | String |
| description | short description of the diagnosis or procedure | String |
| synonyms | List of alternative names and identifiers | List of Strings |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Clinical_variable | HAS_PARENT | Clinical_variable | N/A |
| Phenotype | MAPS_TO | Clinical_variable | N/A |
| Disease | MAPS_TO | Clinical_variable | N/A |
| Biological_sample | HAS_DISEASE | Clinical_variable | genetically_confirmed : String (possible values: True/False) |
| Clinical_variable | REPLACED_BY | Clinical_variable | N/A | 


#### Project

Top level node for the study.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| acronym | Name of the project | String |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Project | HAS_ENROLLED | Subject | N/A |


#### Subject

Node representation of the individual study participant.

##### Properties
| Property | Description | Datatype |
| --- | --- | --- |
| external_id | external id fof the subject for cross-referencing <br> The external_id is "STUDY" for a subjects in the graph | String |
| subjectid | unique number assigned to each subject | int |


##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Project | HAS_ENROLLED | Subject | N/A |
| Biological_sample | BELONGS_TO_SUBJECT | Subject | N/A |


#### Biological_sample

Representation of all samples taken from the subject at a distinct point in time.

##### Properties

| Property | Description | Datatype |
| --- | --- | --- |
| external_id | external id fof the subject for cross-referencing <br> The external_id is "STUDY" for a subjects in the graph | String |
| subjectid | unique number assigned to each subject | int |

##### Relations

| Start | Type | End | Attributes |
| --- | --- | --- | --- |
| Biological_sample | BELONGS_TO_SUBJECT | Subject | N/A |
| Biological_sample | HAS_PHENOTYPE | Phenotype | N/A |
| Biological_sample | HAS_QUANTIFIED_PROTEIN | Protein | source : String (possible values: urine/blood) <br> score : float (log10 quantification) |
| Biological_sample | HAS_DAMAGE | Gene | cadd : String |
| Biological_sample | HAS_DISEASE | Disease | genetically_confirmed : String (possible values: True/False) |
| Biological_sample | HAS_DISEASE | Clinical_variable | genetically_confirmed : String (possible values: True/False) |

