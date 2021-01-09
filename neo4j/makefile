DIR=$(HOME)/github/kg-reactome/neo4j
PROJECT_DIR=$(HOME)/github/kg-geneOntology/neo4j

default:

start:
	docker run \
          -d \
          --name tiny \
          --publish=7575:7474 \
          --publish=7788:7687 \
          --volume=$(DIR)/data:/data \
          --volume=$(DIR)/logs:/logs \
          --volume=$(DIR)/conf:/conf \
          --volume=$(PROJECT_DIR)/import:/import \
          --env=NEO4J_AUTH=none \
          --env=NEO4J_dbms.allow_upgrade=true \
          --env=NEO4J_dbms_memory_heap_max__size=1G \
          neo4j:latest

#          --env=NEO4J_dbms.default_database=fromv3 \
#          --env='NEO4JLABS_PLUGINS=["apoc", "graph-data-science"]' \

bash:
	docker exec -it tiny bash 


stop:
	- docker stop tiny; docker rm tiny;

log:
	docker logs tiny


lsof:
	- lsof -i -n -P | grep 7474
