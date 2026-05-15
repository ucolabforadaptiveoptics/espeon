# %%
# 67.22 micron MMF cladding
# 15.64 micron radius MMF core (eyeballed off image)
# 6.5 micron SMF core
# 125 micron SMF core-to-core
# 10x taper factor (67.22 micron MMF cladding to 633.6 micron SMF jacket)
# 10cm taper distance (arbitrary)

# n_core = 1.45213 and n_clad = 1.44692 for SMF28 https://iopscience.iop.org/article/10.1088/1742-6596/1655/1/012160/pdf
# sqrt(1.44692^2 - n_jacket^2) = 0.13 (design NA from Caleb) -> n_jacket = 1.4410

import numpy as np
from espeon.design import save_lantern_design, ngon_pattern
port_positions_asymmetric = np.load("/Users/adityasengupta/projects/umbreon/data/florida_pl_microscope/florida_asymmetric_centroids.npy")
port_positions_symmetric = ngon_pattern(6, 3, 6.256)

"""wavelengths_nm = np.array([1507.56260406, 1508.58426467, 1509.60401389, 1510.62181916,
       1511.63769679, 1512.65166293, 1513.66368522, 1514.67377986,
       1515.68196292, 1516.68820222, 1517.69251387, 1518.69491386,
       1519.69537018, 1520.69389884, 1521.69051575, 1522.68518908,
       1523.67793476, 1524.66876859, 1525.65765894, 1526.64462163,
       1527.62967238, 1528.61277974, 1529.59395945, 1530.57322712,
       1531.5505515 , 1532.52594822, 1533.49943281, 1534.4709742 ,
       1535.44058794, 1536.40828946, 1537.37404786, 1538.33787861,
       1539.29979705, 1540.25977246, 1541.21782023, 1542.17395559,
       1543.12814802, 1544.0804128 , 1545.03074993, 1545.97917453,
       1546.92565632, 1547.87021046, 1548.81285199, 1549.75355079,
       1550.69232195, 1551.62918039, 1552.56409622, 1553.49708439,
       1554.42815975, 1555.35729259, 1556.28449777, 1557.20979006,
       1558.13313991, 1559.05456211, 1559.97407132, 1560.89163818,
       1561.8072774 , 1562.72100353, 1563.6327874 , 1564.54264363,
       1565.45058668, 1566.35658758, 1567.26066082, 1568.16282079,
       1569.0630387 , 1569.96132896, 1570.85770585, 1571.75214078,
       1572.64464805, 1573.53524186, 1574.4238938 , 1575.31061808,
       1576.19542882, 1577.07829777, 1577.95923907, 1578.83826673,
       1579.7153527 , 1580.59051101])"""

for (tag, pos) in zip(["symmetric", "asymmetric"], [port_positions_asymmetric, port_positions_symmetric]):
    save_lantern_design(
        "florida_" + tag, 
        port_positions=pos, 
        core_radius_um=6.5/2, cladding_radius_um=15.64, z_extent_um=10_000, scale=633.6/(15.64*2),
        n_clad=1.44692, n_core=1.45213, n_jacket=1.441, 
        wavelengths_um=wavelengths_nm/1e3,
        do_plot=True,
        force_overwrite=True,
        simulation_params = {
            "name" : "lightbeam",
            "mesh_extent_um" : 700,
            "mesh_spacing_um" : 0.25,
            "dz_um" : 50,
            "PML" : 8
        }
    )


# %%
