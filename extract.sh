cd $(cat data_directory.txt)

# extract routable data as .sql files using osm2po
java -jar osm2po-core-5.5.2-signed.jar cmd=tjspg prefix=$1 $1-latest.osm.pbf 

# create a new database and activate the postgis extension
sudo service postgresql start
sudo su - postgres -c "dropdb $1"
sudo su - postgres -c "createdb $1 -O EVCC"
sudo su - postgres -c "psql -d $1 -c 'create extension postgis;'"

# load the vertices and edges
cd $1
psql -d $1 -c "\i $1_2po_vertex.sql;"
psql -d $1 -c "\i $1_2po_4pgr.sql;"

# extract data in .csv format
vertex_columns="id, clazz, osm_id, osm_name, ref_count, restrictions, geom_vertex"
edge_columns="osm_id, osm_source_id, osm_target_id, km, kmh, cost, reverse_cost, x1, y1, x2, y2"
psql -d $1 -c "\copy (select $vertex_columns from $1_2po_vertex) to stdout with (format csv, header)" > vertices.csv
psql -d $1 -c "\copy (select $edge_columns from $1_2po_4pgr) to stdout with (format csv, header)" > edges.csv