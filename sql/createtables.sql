CREATE SCHEMA postgis;
GRANT USAGE ON schema postgis to public;
CREATE EXTENSION postgis SCHEMA postgis;
ALTER DATABASE geoprog SET search_path=public,postgis,contrib;

CREATE SCHEMA main

CREATE TABLE main.stationdata 
(
	gid serial primary key,
	geom geometry(point, 4326) NOT NULL,
	intensidadeVentoKM numeric (3, 1) NOT NULL,
	temperatura numeric (3, 1) NOT NULL,
	idEstacao integer NOT NULL,
	pressao numeric (5, 1) NOT NULL,
	humidade numeric (4, 1) NOT NULL,
	localEstacao varchar(80) NOT NULL,
	precAcumulada numeric (3, 1) NOT NULL,
	idDireccVento integer NOT NULL,
	radiacao numeric (5, 1) NOT NULL,
	time timestamp NOT NULL,
	intensidadeVento numeric (3, 1) NOT NULL,
	descDirVento character(10) NOT NULL,
	UNIQUE (idEstacao, time) -- Prevents duplication of data
);

CREATE INDEX stationdata_geom_idx
	ON main.stationdata USING gist(geom);

CREATE INDEX stationdata_idEstacao_idx
	ON main.stationdata(idEstacao);