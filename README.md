## Setup

### Compilation
This project uses cmake. Make sure you have cmake installed by checking:
```bash
cmake --version
```
If you don't have cmake installed, you can do so by running:

```bash
sudo apt install cmake
```
To compile the project run:
```bash
cmake . -B build
cmake --build build --target Routing --config Release
```

You can download our preprocessed graph data as a folder called data from our [GoogleDrive](https://drive.google.com/drive/folders/1KUM74qZA9uC3vnbOFeFcUWTnsZufuuS1?usp=sharing) to an arbitrary location of your choice.

Make sure to copy the path to your data folder into **data_directory.txt**.
The path must be relative to the root folder of the repository.

## Basic Usage
Program options can be displayed with *-h (--help)*.
| option | arg | description |
|---|---|---|
| -h, --help |  | Print help |
| -g, --graph | string | Input graph (e.g. germany, ...) |
| -r, --runs | number | Number of test queries (default: 0) |
| -p, --partition_size | number | Partition size (default: 128) |
| -t, --threads | number | Number of parallel threads (default: 4) |
| -v, --validate |  | Validate query results against dijkstra |
| -f, --force |  | Force flag precomputation (overwrite existing flags) |

## Graph creation
### Java

A java runtime environment is required.

```bash
sudo apt install default-jre
```

### osm2po

Go to [osm2po.de](https://osm2po.de) and download osm2po version 5.5.2. located in the ZIP file in the download section of the webpage. Extract the whole ZIP content into the data folder.

### PostgreSQL and PostGIS

Install PostgreSQL and PostGIS! To achieve that just copy the following commands into your terminal:

```bash
sudo apt -y install gnupg2; 
wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -; 
echo "deb http://apt.postgresql.org/pub/repos/apt/ `lsb_release -cs`-pgdg main" | sudo tee /etc/apt/sources.list.d/pgdg.list;
```

```bash
sudo apt -y install postgresql-12 postgresql-client-12;
```

```bash
sudo apt install postgis postgresql-12-postgis-3;
```

```bash
sudo apt-get install postgresql-12-postgis-3-scripts 
```

Now create a user called "EVCC". (Another script uses this name.)

```bash
sudo su - postgres -c "createuser username"
```

Also execute:

```bash
sudo su - postgres -c "createuser your_username"
```

Finally replace the config file in the data folder by RoutingPipeline/osm2po.config.

## Graph Data Input

You can download OpenStreeMap data from the private company [Geofabrik](http://download.geofabrik.de/) from Karlsruhe.
Put it into the data folder.

Now extract the data via
```sh extract.sh region```
where region-latest.osm.pbf is the file you downloaded.
