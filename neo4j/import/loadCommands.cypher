LOAD CSV with HEADERS from 'file:///ensgSymbols.tsv' AS line FIELDTERMINATOR '\t'
     CREATE (:Entity {id:   line.GENEID,
                      name: line.SYMBOL,
                      type: line.type});

LOAD CSV with HEADERS from 'file:///uniprotSymbols.tsv' AS line FIELDTERMINATOR '\t'
     CREATE (:Entity {id:   line.UNIPROTID,
                      name: line.SYMBOL,
                      type: line.type});

LOAD CSV with HEADERS from 'file:///goidTermsOntology.tsv' AS line FIELDTERMINATOR '\t'
     CREATE (:Entity {id:   line.GOID,
                      name: line.TERM,
		      ontology: line.ONTOLOGY,
                      type: line.type});

DROP INDEX ON :Entity(id);
CREATE INDEX ON :Entity(id);



LOAD CSV WITH HEADERS FROM  "file:///ensgGO.tsv" AS line FIELDTERMINATOR '\t'
      MATCH (source:Entity {id: line.GO}), (target:Entity {id: line.GENEID})
      CREATE (source)-[:interacts {type: line.interaction,
      	                           ontology: line.ONTOLOGY,
				   evidenece: line.EVIDENCE}]->(target);
      

