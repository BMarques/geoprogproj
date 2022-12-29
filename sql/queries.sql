-- Get all last measured temps
SELECT geom, temperatura::float, localestacao
FROM main.stationdata
WHERE time = 
(
    SELECT max(time)
    FROM main.stationdata
)
AND temperatura <> -99.0