* Encoding: UTF-8.
 
****Compute composites of CBCL raw items for each wave based on methods of Michelini et al. (2019)****

COMPUTE comp1_1=SUM(cbcl_q57_p1,cbcl_q97_p1).
EXECUTE.
COMPUTE comp1_2=SUM(cbcl_q57_p2,cbcl_q97_p2).
EXECUTE.
COMPUTE comp1_3=SUM(cbcl_q57_p3,cbcl_q97_p3).
EXECUTE.

COMPUTE comp2_1=SUM(cbcl_q22_p1,cbcl_q23_p1,cbcl_q28_p1).
EXECUTE.
COMPUTE comp2_2=SUM(cbcl_q22_p2,cbcl_q23_p2,cbcl_q28_p2).
EXECUTE.
COMPUTE comp2_3=SUM(cbcl_q22_p3,cbcl_q23_p3,cbcl_q28_p3).
EXECUTE.

COMPUTE comp3_1=SUM(cbcl_q20_p1,cbcl_q21_p1,cbcl_q106_p1).
EXECUTE.
COMPUTE comp3_2=SUM(cbcl_q20_p2,cbcl_q21_p2,cbcl_q106_p2).
EXECUTE.
COMPUTE comp3_3=SUM(cbcl_q20_p3,cbcl_q21_p3,cbcl_q106_p3).
EXECUTE.

COMPUTE comp4_1=SUM(cbcl_q81_p1,cbcl_q82_p1).
EXECUTE.
COMPUTE comp4_2=SUM(cbcl_q81_p2,cbcl_q82_p2).
EXECUTE.
COMPUTE comp4_3=SUM(cbcl_q81_p3,cbcl_q82_p3).
EXECUTE.

COMPUTE comp5_1=SUM(cbcl_q25_p1,cbcl_q48_p1).
EXECUTE.
COMPUTE comp5_2=SUM(cbcl_q25_p2,cbcl_q48_p2).
EXECUTE.
COMPUTE comp5_3=SUM(cbcl_q25_p3,cbcl_q48_p3).
EXECUTE.

COMPUTE comp6_1=SUM(cbcl_q08_p1,cbcl_q78_p1,cbcl_q10_p1).
EXECUTE.
COMPUTE comp6_2=SUM(cbcl_q08_p2,cbcl_q78_p2,cbcl_q10_p2).
EXECUTE.
COMPUTE comp6_3=SUM(cbcl_q08_p3,cbcl_q78_p3,cbcl_q10_p3).
EXECUTE. 
 
****Create unstandardized factor scores for EXT, INT, ND, SOM, DET, and p from higher-order model and EXT, INT, ND, SOM, DET from correlated factors model for baseline wave 1, wave 2, & wave 3*******
****CBCL item-level variables in this syntax were taken from raw ABCD file "abcd_cbcl01.txt"*****

****Wave 1****

****Multiply CBCL items by unstandardized weight from higher-order model****

****EXT factor from higher-order model

COMPUTE comp1_weight1=comp1_1 * 0.441.
EXECUTE.
COMPUTE q16_weight1=cbcl_q16_p1 * 0.411.
EXECUTE.
COMPUTE comp2_weight1=comp2_1 * 0.456.
EXECUTE.
COMPUTE q37_weight1=cbcl_q37_p1 * 0.412.
EXECUTE.
COMPUTE q95_weight1=cbcl_q95_p1 * 0.444.
EXECUTE.
COMPUTE q3_weight1=cbcl_q03_p1 * 0.438.
EXECUTE.
COMPUTE comp3_weight1=comp3_1 * 0.438.
EXECUTE.
COMPUTE q68_weight1=cbcl_q68_p1 * 0.422.
EXECUTE.
COMPUTE q26_weight1=cbcl_q26_p1 * 0.401.
EXECUTE.
COMPUTE q90_weight1=cbcl_q90_p1 * 0.361.
EXECUTE.
COMPUTE q94_weight1=cbcl_q94_p1 * 0.379.
EXECUTE.
COMPUTE comp4_weight1=comp4_1 * 0.391.
EXECUTE.
COMPUTE q86_weight1=cbcl_q86_p1 * 0.440.
EXECUTE.
COMPUTE q43_weight1=cbcl_q43_p1 * 0.396.
EXECUTE.
COMPUTE q67_weight1=cbcl_q67_p1 * 0.410.
EXECUTE.
COMPUTE q87_weight1=cbcl_q87_p1 * 0.453.
EXECUTE.
COMPUTE q27_weight1=cbcl_q27_p1 * 0.395.
EXECUTE.
COMPUTE comp5_weight1=comp5_1 * 0.429.
EXECUTE.
COMPUTE q89_weight1=cbcl_q89_p1 * 0.420.
EXECUTE.
COMPUTE q19_weight1=cbcl_q19_p1 * 0.419.
EXECUTE.
COMPUTE q39_weight1=cbcl_q39_p1 * 0.314.
EXECUTE.
COMPUTE q34_weight1=cbcl_q34_p1 * 0.410.
EXECUTE.
COMPUTE q88_weight1=cbcl_q88_p1 * 0.417.
EXECUTE.
COMPUTE q7_weight1=cbcl_q07_p1 * 0.304.
EXECUTE.
COMPUTE q109_weight1=cbcl_q109_p1 * 0.361.
EXECUTE.

****Compute unstandardized factor scores for EXT from higher-order model

COMPUTE EXT_fs_ho1=SUM(comp1_weight1,q16_weight1,comp2_weight1,q37_weight1,q95_weight1,comp3_weight1,
    q68_weight1,q26_weight1,q90_weight1,q94_weight1,comp4_weight1,q86_weight1,q43_weight1,q67_weight1,
    q87_weight1,q27_weight1,comp5_weight1,q89_weight1,q19_weight1,q39_weight1,q34_weight1,q88_weight1,
    q109_weight1,q3_weight1,q7_weight1).
EXECUTE.

****INT factor from higher-order model

COMPUTE q50_weight1=cbcl_q50_p1 * 0.469.
EXECUTE.
COMPUTE q112_weight1=cbcl_q112_p1 * 0.441.
EXECUTE.
COMPUTE q32_weight1=cbcl_q32_p1 * 0.299.
EXECUTE.
COMPUTE q52_weight1=cbcl_q52_p1 * 0.425.
EXECUTE.
COMPUTE q45_weight1=cbcl_q45_p1 * 0.492.
EXECUTE.
COMPUTE q31_weight1=cbcl_q31_p1 * 0.403.
EXECUTE.
COMPUTE q35_weight1=cbcl_q35_p1 * 0.478.
EXECUTE.
COMPUTE q30_weight1=cbcl_q30_p1 * 0.410.
EXECUTE.
COMPUTE q29_weight1=cbcl_q29_p1 * 0.335.
EXECUTE.
COMPUTE q12_weight1=cbcl_q12_p1 * 0.455.
EXECUTE.

****Compute unstandardized factor scores for INT from higher-order model

COMPUTE INT_fs_ho1=SUM(q50_weight1,q112_weight1,q32_weight1,q52_weight1,q45_weight1,q31_weight1,q35_weight1,
    q30_weight1,q29_weight1,q12_weight1).
EXECUTE.

****ND factor from higher-order model

COMPUTE comp6_weight1=comp6_1 * 0.351.
EXECUTE.
COMPUTE q17_weight1=cbcl_q17_p1 * 0.266.
EXECUTE.
COMPUTE q13_weight1=cbcl_q13_p1 * 0.310.
EXECUTE.
COMPUTE q62_weight1=cbcl_q62_p1 * 0.292.
EXECUTE.
COMPUTE q46_weight1=cbcl_q46_p1 * 0.280.
EXECUTE.
COMPUTE q4_weight1=cbcl_q04_p1 * 0.336.
EXECUTE.
COMPUTE q9_weight1=cbcl_q09_p1 * 0.335.
EXECUTE.
COMPUTE q61_weight1=cbcl_q61_p1 * 0.308.
EXECUTE.
COMPUTE q66_weight1=cbcl_q66_p1 * 0.325.
EXECUTE.
COMPUTE q1_weight1=cbcl_q01_p1 * 0.280.
EXECUTE.
COMPUTE q36_weight1=cbcl_q36_p1 * 0.239.
EXECUTE.
COMPUTE q64_weight1=cbcl_q64_p1 * 0.246.
EXECUTE.

****Compute unstandardized factor scores for ND from higher-order model

COMPUTE ND_fs_ho1=SUM(comp6_weight1,q17_weight1,q13_weight1,q62_weight1,q46_weight1,q4_weight1,q9_weight1,
    q61_weight1,q66_weight1,q1_weight1,q36_weight1,q64_weight1).
EXECUTE.

****SOM factor from higher-order model

COMPUTE q56c_weight1=cbcl_q56c_p1 * 0.694.
EXECUTE.
COMPUTE q56f_weight1=cbcl_q56f_p1 * 0.623.
EXECUTE.
COMPUTE q56g_weight1=cbcl_q56g_p1 * 0.501.
EXECUTE.
COMPUTE q56b_weight1=cbcl_q56b_p1 * 0.527.
EXECUTE.
COMPUTE q56a_weight1=cbcl_q56a_p1 * 0.532.
EXECUTE.
COMPUTE q51_weight1=cbcl_q51_p1 * 0.641.
EXECUTE.
COMPUTE q56d_weight1=cbcl_q56d_p1 * 0.439.
EXECUTE.
COMPUTE q56e_weight1=cbcl_q56e_p1 * 0.403.
EXECUTE.

****Compute unstandardized factor scores for SOM from higher-order model

COMPUTE SOM_fs_ho1=SUM(q56c_weight1,q56f_weight1,q56g_weight1,q56b_weight1,q56a_weight1,q51_weight1,
    q56d_weight1,q56e_weight1).
EXECUTE.

****DET factor from higher-order model

COMPUTE q111_weight1=cbcl_q111_p1 * 0.566.
EXECUTE.
COMPUTE q42_weight1=cbcl_q42_p1 * 0.446.
EXECUTE.
COMPUTE q75_weight1=cbcl_q75_p1 * 0.366.
EXECUTE.
COMPUTE q65_weight1=cbcl_q65_p1 * 0.509.
EXECUTE.
COMPUTE q102_weight1=cbcl_q102_p1 * 0.491.
EXECUTE.

****Compute unstandardized factor scores for DET from higher-order model

COMPUTE DET_fs_ho1=SUM(q111_weight1,q42_weight1,q75_weight1,q65_weight1,q102_weight1).
EXECUTE.

****p factor from higher-order model

COMPUTE EXT_fs_ho_weight1=EXT_fs_ho3 * 1.465.
EXECUTE.
COMPUTE INT_fs_ho_weight1=INT_fs_ho3 * 1.357.
EXECUTE.
COMPUTE ND_fs_ho_weight1=ND_fs_ho3 * 2.069.
EXECUTE.
COMPUTE SOM_fs_ho_weight1=SOM_fs_ho3 * 0.691.
EXECUTE.
COMPUTE DET_fs_ho_weight1=DET_fs_ho3 * 1.169.
EXECUTE.

****Compute unstandardized factor scores for p from higher-order model

COMPUTE p_fs_ho1=SUM(EXT_fs_ho_weight1,INT_fs_ho_weight1,ND_fs_ho_weight1,SOM_fs_ho_weight1,
    DET_fs_ho_weight1).
EXECUTE.


****Multiply CBCL items by unstandardized weight1 from correlated factors model****

****EXT factor from correlated factors model

COMPUTE comp1_corr_weight1=comp1_1 * 0.782.
EXECUTE.
COMPUTE q16_corr_weight1=cbcl_q16_p1 * 0.729.
EXECUTE.
COMPUTE comp2_corr_weight1=comp2_1 * 0.809.
EXECUTE.
COMPUTE q37_corr_weight1=cbcl_q37_p1 * 0.731.
EXECUTE.
COMPUTE q95_corr_weight1=cbcl_q95_p1 * 0.787.
EXECUTE.
COMPUTE q3_corr_weight1=cbcl_q03_p1 * 0.776.
EXECUTE.
COMPUTE comp3_corr_weight1=comp3_1 * 0.778.
EXECUTE.
COMPUTE q68_corr_weight1=cbcl_q68_p1 * 0.749.
EXECUTE.
COMPUTE q26_corr_weight1=cbcl_q26_p1 * 0.711.
EXECUTE.
COMPUTE q90_corr_weight1=cbcl_q90_p1 * 0.641.
EXECUTE.
COMPUTE q94_corr_weight1=cbcl_q94_p1 * 0.672.
EXECUTE.
COMPUTE comp4_corr_weight1=comp4_1 * 0.694.
EXECUTE.
COMPUTE q86_corr_weight1=cbcl_q86_p1 * 0.781.
EXECUTE.
COMPUTE q43_corr_weight1=cbcl_q43_p1 * 0.703.
EXECUTE.
COMPUTE q67_corr_weight1=cbcl_q67_p1 * 0.728.
EXECUTE.
COMPUTE q87_corr_weight1=cbcl_q87_p1 * 0.803.
EXECUTE.
COMPUTE q27_corr_weight1=cbcl_q27_p1 * 0.701.
EXECUTE.
COMPUTE comp5_corr_weight1=comp5_1 * 0.760.
EXECUTE.
COMPUTE q89_corr_weight1=cbcl_q89_p1 * 0.744.
EXECUTE.
COMPUTE q19_corr_weight1=cbcl_q19_p1 * 0.744.
EXECUTE.
COMPUTE q39_corr_weight1=cbcl_q39_p1 * 0.558.
EXECUTE.
COMPUTE q34_corr_weight1=cbcl_q34_p1 * 0.726.
EXECUTE.
COMPUTE q88_corr_weight1=cbcl_q88_p1 * 0.739.
EXECUTE.
COMPUTE q7_corr_weight1=cbcl_q07_p1 * 0.539.
EXECUTE.
COMPUTE q109_corr_weight1=cbcl_q109_p1 * 0.640.
EXECUTE.

****Compute unstandardized factor scores for EXT from correlated factors model

COMPUTE EXT_fs_corr1=SUM(comp1_corr_weight1,q16_corr_weight1,comp2_corr_weight1,q37_corr_weight1,
    q95_corr_weight1,q3_corr_weight1,comp3_corr_weight1,q68_corr_weight1,q26_corr_weight1,q90_corr_weight1,
    q94_corr_weight1,comp4_corr_weight1,q86_corr_weight1,q43_corr_weight1,q67_corr_weight1,q87_corr_weight1,
    q27_corr_weight1,comp5_corr_weight1,q89_corr_weight1,q19_corr_weight1,q39_corr_weight1,q34_corr_weight1,
    q88_corr_weight1,q7_corr_weight1,q109_corr_weight1).
EXECUTE.

****INT factor from correlated factors model

COMPUTE q50_corr_weight1=cbcl_q50_p1 * 0.791.
EXECUTE.
COMPUTE q112_corr_weight1=cbcl_q112_p1 * 0.745.
EXECUTE.
COMPUTE q32_corr_weight1=cbcl_q32_p1 * 0.505.
EXECUTE.
COMPUTE q52_corr_weight1=cbcl_q52_p1 * 0.717.
EXECUTE.
COMPUTE q45_corr_weight1=cbcl_q45_p1 * 0.828.
EXECUTE.
COMPUTE q31_corr_weight1=cbcl_q31_p1 * 0.678.
EXECUTE.
COMPUTE q35_corr_weight1=cbcl_q35_p1 * 0.804.
EXECUTE.
COMPUTE q30_corr_weight1=cbcl_q30_p1 * 0.692.
EXECUTE.
COMPUTE q29_corr_weight1=cbcl_q29_p1 * 0.564.
EXECUTE.
COMPUTE q12_corr_weight1=cbcl_q12_p1 * 0.765.
EXECUTE.

****Compute unstandardized factor scores for INT from correlated factors model

COMPUTE INT_fs_corr1=SUM(q50_corr_weight1,q112_corr_weight1,q32_corr_weight1,q52_corr_weight1,
    q45_corr_weight1,q31_corr_weight1,q35_corr_weight1,q30_corr_weight1,q29_corr_weight1,q12_corr_weight1).
EXECUTE.

****ND factor from correlated factors model

COMPUTE comp6_corr_weight1=comp6_1 * 0.808.
EXECUTE.
COMPUTE q17_corr_weight1=cbcl_q17_p1 * 0.612.
EXECUTE.
COMPUTE q13_corr_weight1=cbcl_q13_p1 * 0.712.
EXECUTE.
COMPUTE q62_corr_weight1=cbcl_q62_p1 * 0.671.
EXECUTE.
COMPUTE q46_corr_weight1=cbcl_q46_p1 * 0.644.
EXECUTE.
COMPUTE q4_corr_weight1=cbcl_q04_p1 * 0.772.
EXECUTE.
COMPUTE q9_corr_weight1=cbcl_q09_p1 * 0.769.
EXECUTE.
COMPUTE q61_corr_weight1=cbcl_q61_p1 * 0.708.
EXECUTE.
COMPUTE q66_corr_weight1=cbcl_q66_p1 * 0.748.
EXECUTE.
COMPUTE q1_corr_weight1=cbcl_q01_p1 * 0.644.
EXECUTE.
COMPUTE q36_corr_weight1=cbcl_q36_p1 * 0.548.
EXECUTE.
COMPUTE q64_corr_weight1=cbcl_q64_p1 * 0.565.
EXECUTE.

****Compute unstandardized factor scores for ND from correlated factors model

COMPUTE ND_fs_corr1=SUM(comp6_corr_weight1,q17_corr_weight1,q13_corr_weight1,q62_corr_weight1,
    q46_corr_weight1,q4_corr_weight1,q9_corr_weight1,q61_corr_weight1,q66_corr_weight1,q1_corr_weight1,
    q36_corr_weight1,q64_corr_weight1).
EXECUTE.

****SOM factor from correlated factors model

COMPUTE q56c_corr_weight1=cbcl_q56c_p1 * 0.844.
EXECUTE.
COMPUTE q56f_corr_weight1=cbcl_q56f_p1 * 0.759.
EXECUTE.
COMPUTE q56g_corr_weight1=cbcl_q56g_p1 * 0.607.
EXECUTE.
COMPUTE q56b_corr_weight1=cbcl_q56b_p1 * 0.641.
EXECUTE.
COMPUTE q56a_corr_weight1=cbcl_q56a_p1 * 0.646.
EXECUTE.
COMPUTE q51_corr_weight1=cbcl_q51_p1 * 0.779.
EXECUTE.
COMPUTE q56d_corr_weight1=cbcl_q56d_p1 * 0.532.
EXECUTE.
COMPUTE q56e_corr_weight1=cbcl_q56e_p1 * 0.489.
EXECUTE.

****Compute unstandardized factor scores for SOM from correlated factors model

COMPUTE SOM_fs_corr1=SUM(q56c_corr_weight1,q56f_corr_weight1,q56g_corr_weight1,q56b_corr_weight1,
    q56a_corr_weight1,q51_corr_weight1,q56d_corr_weight1,q56e_corr_weight1).
EXECUTE.

****DET factor from correlated factors model

COMPUTE q111_corr_weight1=cbcl_q111_p1 * 0.869.
EXECUTE.
COMPUTE q42_corr_weight1=cbcl_q42_p1 * 0.686.
EXECUTE.
COMPUTE q75_corr_weight1=cbcl_q75_p1 * 0.567.
EXECUTE.
COMPUTE q65_corr_weight1=cbcl_q65_p1 * 0.781.
EXECUTE.
COMPUTE q102_corr_weight1=cbcl_q102_p1 * 0.754.
EXECUTE.

****Compute unstandardized factor scores for DET from correlated factors model

COMPUTE DET_fs_corr1=SUM(q111_corr_weight1,q42_corr_weight1,q75_corr_weight1,q65_corr_weight1,
    q102_corr_weight1).
EXECUTE.

****Multiply CBCL items by unstandardized weight1 from the bifactor model****

****EXT factor from bifactor model

COMPUTE comp1_bi_ext_weight1=comp1_1 * 0.628.
EXECUTE.
COMPUTE q16_bi_ext_weight1=cbcl_q16_p1 * 0.628.
EXECUTE.
COMPUTE comp2_bi_ext_weight1=comp2_1 * 0.558.
EXECUTE.
COMPUTE q37_bi_ext_weight1=cbcl_q37_p1 * 0.577.
EXECUTE.
COMPUTE q95_bi_ext_weight1=cbcl_q95_p1 * 0.473.
EXECUTE.
COMPUTE q3_bi_ext_weight1=cbcl_q03_p1 * 0.489.
EXECUTE.
COMPUTE comp3_bi_ext_weight1=comp3_1 * 0.484.
EXECUTE.
COMPUTE q68_bi_ext_weight1=cbcl_q68_p1 * 0.449.
EXECUTE.
COMPUTE q26_bi_ext_weight1=cbcl_q26_p1 * 0.502.
EXECUTE.
COMPUTE q90_bi_ext_weight1=cbcl_q90_p1 * 0.481.
EXECUTE.
COMPUTE q94_bi_ext_weight1=cbcl_q94_p1 * 0.493.
EXECUTE.
COMPUTE comp4_bi_ext_weight1=comp4_1 * 0.548.
EXECUTE.
COMPUTE q86_bi_ext_weight1=cbcl_q86_p1 * 0.364.
EXECUTE.
COMPUTE q43_bi_ext_weight1=cbcl_q43_p1 * 0.514.
EXECUTE.
COMPUTE q67_bi_ext_weight1=cbcl_q67_p1 * 0.325.
EXECUTE.
COMPUTE q87_bi_ext_weight1=cbcl_q87_p1 * 0.253.
EXECUTE.
COMPUTE q27_bi_ext_weight1=cbcl_q27_p1 * 0.312.
EXECUTE.
COMPUTE comp5_bi_ext_weight1=comp5_1 * 0.257.
EXECUTE.
COMPUTE q89_bi_ext_weight1=cbcl_q89_p1 * 0.267.
EXECUTE.
COMPUTE q19_bi_ext_weight1=cbcl_q19_p1 * 0.273.
EXECUTE.
COMPUTE q39_bi_ext_weight1=cbcl_q39_p1 * 0.405.
EXECUTE.
COMPUTE q34_bi_ext_weight1=cbcl_q34_p1 * 0.150.
EXECUTE.
COMPUTE q88_bi_ext_weight1=cbcl_q88_p1 * 0.156.
EXECUTE.
COMPUTE q7_bi_ext_weight1=cbcl_q07_p1 * 0.370.
EXECUTE.
COMPUTE q109_bi_ext_weight1=cbcl_q109_p1 * 0.189.
EXECUTE.

****Compute unstandardized factor scores for EXT from  bifactor model

COMPUTE EXT_fs_bi1=SUM(comp1_bi_ext_weight1,q16_bi_ext_weight1,comp2_bi_ext_weight1,q37_bi_ext_weight1,
    q95_bi_ext_weight1,q3_bi_ext_weight1,comp3_bi_ext_weight1,q68_bi_ext_weight1,q26_bi_ext_weight1,q90_bi_ext_weight1,
    q94_bi_ext_weight1,comp4_bi_ext_weight1,q86_bi_ext_weight1,q43_bi_ext_weight1,q67_bi_ext_weight1,q87_bi_ext_weight1,
    q27_bi_ext_weight1,comp5_bi_ext_weight1,q89_bi_ext_weight1,q19_bi_ext_weight1,q39_bi_ext_weight1,q34_bi_ext_weight1,
    q88_bi_ext_weight1,q7_bi_ext_weight1,q109_bi_ext_weight1).
EXECUTE.

****INT factor from bifactor model

COMPUTE q50_bi_int_weight1=cbcl_q50_p1 * 0.599.
EXECUTE.
COMPUTE q112_bi_int_weight1=cbcl_q112_p1 * 0.555.
EXECUTE.
COMPUTE q32_bi_int_weight1=cbcl_q32_p1 * 0.538.
EXECUTE.
COMPUTE q52_bi_int_weight1=cbcl_q52_p1 * 0.555.
EXECUTE.
COMPUTE q45_bi_int_weight1=cbcl_q45_p1 * 0.395.
EXECUTE.
COMPUTE q31_bi_int_weight1=cbcl_q31_p1 * 0.476.
EXECUTE.
COMPUTE q35_bi_int_weight1=cbcl_q35_p1 * 0.301.
EXECUTE.
COMPUTE q30_bi_int_weight1=cbcl_q30_p1 * 0.298.
EXECUTE.
COMPUTE q29_bi_int_weight1=cbcl_q29_p1 * 0.333.
EXECUTE.
COMPUTE q12_bi_int_weight1=cbcl_q12_p1 * 0.129.
EXECUTE.

****Compute unstandardized factor scores for INT from bifactor model

COMPUTE INT_fs_bi1=SUM(q50_bi_int_weight1,q112_bi_int_weight1,q32_bi_int_weight1,q52_bi_int_weight1,
    q45_bi_int_weight1,q31_bi_int_weight1,q35_bi_int_weight1,q30_bi_int_weight1,q29_bi_int_weight1,q12_bi_int_weight1).
EXECUTE.

****ND factor from bifactor model

COMPUTE comp6_bi_nd_weight1=comp6_1 * 0.489.
EXECUTE.
COMPUTE q17_bi_nd_weight1=cbcl_q17_p1 * 0.413.
EXECUTE.
COMPUTE q13_bi_nd_weight1=cbcl_q13_p1 * 0.383.
EXECUTE.
COMPUTE q62_bi_nd_weight1=cbcl_q62_p1 * 0.417.
EXECUTE.
COMPUTE q46_bi_nd_weight1=cbcl_q46_p1 * 0.144.
EXECUTE.
COMPUTE q4_bi_nd_weight1=cbcl_q04_p1 * 0.365.
EXECUTE.
COMPUTE q9_bi_nd_weight1=cbcl_q09_p1 * 0.129.
EXECUTE.
COMPUTE q61_bi_nd_weight1=cbcl_q61_p1 * 0.357.
EXECUTE.
COMPUTE q66_bi_nd_weight1=cbcl_q66_p1 * 0.137.
EXECUTE.
COMPUTE q1_bi_nd_weight1=cbcl_q01_p1 * 0.251.
EXECUTE.
COMPUTE q36_bi_nd_weight1=cbcl_q36_p1 * 0.293.
EXECUTE.
COMPUTE q64_bi_nd_weight1=cbcl_q64_p1 * 0.149.
EXECUTE.

****Compute unstandardized factor scores for ND from bifactor model

COMPUTE ND_fs_bi1=SUM(comp6_bi_nd_weight1,q17_bi_nd_weight1,q13_bi_nd_weight1,q62_bi_nd_weight1,
    q46_bi_nd_weight1,q4_bi_nd_weight1,q9_bi_nd_weight1,q61_bi_nd_weight1,q66_bi_nd_weight1,q1_bi_nd_weight1,
    q36_bi_nd_weight1,q64_bi_nd_weight1).
EXECUTE.

****SOM factor from bifactor model

COMPUTE q56c_bi_som_weight1=cbcl_q56c_p1 * 0.776.
EXECUTE.
COMPUTE q56f_bi_som_weight1=cbcl_q56f_p1 * 0.713.
EXECUTE.
COMPUTE q56g_bi_som_weight1=cbcl_q56g_p1 * 0.614.
EXECUTE.
COMPUTE q56b_bi_som_weight1=cbcl_q56b_p1 * 0.533.
EXECUTE.
COMPUTE q56a_bi_som_weight1=cbcl_q56a_p1 * 0.438.
EXECUTE.
COMPUTE q51_bi_som_weight1=cbcl_q51_p1 * 0.443.
EXECUTE.
COMPUTE q56d_bi_som_weight1=cbcl_q56d_p1 * 0.274.
EXECUTE.
COMPUTE q56e_bi_som_weight1=cbcl_q56e_p1 * 0.252.
EXECUTE.

****Compute unstandardized factor scores for SOM from bifactor model

COMPUTE SOM_fs_bi1=SUM(q56c_bi_som_weight1,q56f_bi_som_weight1,q56g_bi_som_weight1,q56b_bi_som_weight1,
    q56a_bi_som_weight1,q51_bi_som_weight1,q56d_bi_som_weight1,q56e_bi_som_weight1).
EXECUTE.

****DET factor from bifactor model

COMPUTE q111_bi_det_weight1=cbcl_q111_p1 * 0.630.
EXECUTE.
COMPUTE q42_bi_det_weight1=cbcl_q42_p1 * 0.531.
EXECUTE.
COMPUTE q75_bi_det_weight1=cbcl_q75_p1 * 0.474.
EXECUTE.
COMPUTE q65_bi_det_weight1=cbcl_q65_p1 * 0.398.
EXECUTE.
COMPUTE q102_bi_det_weight1=cbcl_q102_p1 * 0.299.
EXECUTE.

****Compute unstandardized factor scores for DET from bifactor model

COMPUTE DET_fs_bi1=SUM(q111_bi_det_weight1,q42_bi_det_weight1,q75_bi_det_weight1,q65_bi_det_weight1,
    q102_bi_det_weight1).
EXECUTE.

****p factor from bifactor model

COMPUTE comp1_bi_p_weight1=comp1_1 * 0.558.
EXECUTE.
COMPUTE q16_bi_p_weight1=cbcl_q16_p1 * 0.503.
EXECUTE.
COMPUTE comp2_bi_p_weight1=comp2_1 * 0.626.
EXECUTE.
COMPUTE q37_bi_p_weight1=cbcl_q37_p1 * 0.527.
EXECUTE.
COMPUTE q95_bi_p_weight1=cbcl_q95_p1 * 0.640.
EXECUTE.
COMPUTE q3_bi_p_weight1=cbcl_q03_p1 * 0.622.
EXECUTE.
COMPUTE comp3_bi_p_weight1=comp3_1 * 0.623.
EXECUTE.
COMPUTE q68_bi_p_weight1=cbcl_q68_p1 * 0.609.
EXECUTE.
COMPUTE q26_bi_p_weight1=cbcl_q26_p1 * 0.543.
EXECUTE.
COMPUTE q90_bi_p_weight1=cbcl_q90_p1 * 0.480.
EXECUTE.
COMPUTE q94_bi_p_weight1=cbcl_q94_p1 * 0.507.
EXECUTE.
COMPUTE comp4_bi_p_weight1=comp4_1 * 0.503.
EXECUTE.
COMPUTE q86_bi_p_weight1=cbcl_q86_p1 * 0.679.
EXECUTE.
COMPUTE q43_bi_p_weight1=cbcl_q43_p1 * 0.532.
EXECUTE.
COMPUTE q67_bi_p_weight1=cbcl_q67_p1 * 0.638.
EXECUTE.
COMPUTE q87_bi_p_weight1=cbcl_q87_p1 * 0.745.
EXECUTE.
COMPUTE q27_bi_p_weight1=cbcl_q27_p1 * 0.614.
EXECUTE.
COMPUTE comp5_bi_p_weight1=comp5_1 * 0.699.
EXECUTE.
COMPUTE q89_bi_p_weight1=cbcl_q89_p1 * 0.679.
EXECUTE.
COMPUTE q19_bi_p_weight1=cbcl_q19_p1 * 0.674.
EXECUTE.
COMPUTE q39_bi_p_weight1=cbcl_q39_p1 * 0.424.
EXECUTE.
COMPUTE q34_bi_p_weight1=cbcl_q34_p1 * 0.701.
EXECUTE.
COMPUTE q88_bi_p_weight1=cbcl_q88_p1 * 0.713.
EXECUTE.
COMPUTE q7_bi_p_weight1=cbcl_q07_p1 * 0.420.
EXECUTE.
COMPUTE q109_bi_p_weight1=cbcl_q109_p1 * 0.596.
EXECUTE.
COMPUTE q50_bi_p_weight1=cbcl_q50_p1 * 0.608.
EXECUTE.
COMPUTE q112_bi_p_weight1=cbcl_q112_p1 * 0.573.
EXECUTE.
COMPUTE q32_bi_p_weight1=cbcl_q32_p1 * 0.359.
EXECUTE.
COMPUTE q52_bi_p_weight1=cbcl_q52_p1 * 0.547.
EXECUTE.
COMPUTE q45_bi_p_weight1=cbcl_q45_p1 * 0.683.
EXECUTE.
COMPUTE q31_bi_p_weight1=cbcl_q31_p1 * 0.531.
EXECUTE.
COMPUTE q35_bi_p_weight1=cbcl_q35_p1 * 0.681.
EXECUTE.
COMPUTE q30_bi_p_weight1=cbcl_q30_p1 * 0.578.
EXECUTE.
COMPUTE q29_bi_p_weight1=cbcl_q29_p1 * 0.454.
EXECUTE.
COMPUTE q12_bi_p_weight1=cbcl_q12_p1 * 0.673.
EXECUTE.
COMPUTE comp6_bi_p_weight1=comp6_1 * 0.705.
EXECUTE.
COMPUTE q17_bi_p_weight1=cbcl_q17_p1 * 0.524.
EXECUTE.
COMPUTE q13_bi_p_weight1=cbcl_q13_p1 * 0.628.
EXECUTE.
COMPUTE q62_bi_p_weight1=cbcl_q62_p1 * 0.581.
EXECUTE.
COMPUTE q46_bi_p_weight1=cbcl_q46_p1 * 0.599.
EXECUTE.
COMPUTE q4_bi_p_weight1=cbcl_q04_p1 * 0.685.
EXECUTE.
COMPUTE q9_bi_p_weight1=cbcl_q09_p1 * 0.721.
EXECUTE.
COMPUTE q61_bi_p_weight1=cbcl_q61_p1 * 0.624.
EXECUTE.
COMPUTE q66_bi_p_weight1=cbcl_q66_p1 * 0.700.
EXECUTE.
COMPUTE q1_bi_p_weight1=cbcl_q01_p1 * 0.582.
EXECUTE.
COMPUTE q36_bi_p_weight1=cbcl_q36_p1 * 0.483.
EXECUTE.
COMPUTE q64_bi_p_weight1=cbcl_q64_p1 * 0.522.
EXECUTE.
COMPUTE q56c_bi_p_weight1=cbcl_q56c_p1 * 0.450.
EXECUTE.
COMPUTE q56f_bi_p_weight1=cbcl_q56f_p1 * 0.404.
EXECUTE.
COMPUTE q56g_bi_p_weight1=cbcl_q56g_p1 * 0.310.
EXECUTE.
COMPUTE q56b_bi_p_weight1=cbcl_q56b_p1 * 0.363.
EXECUTE.
COMPUTE q56a_bi_p_weight1=cbcl_q56a_p1 * 0.390.
EXECUTE.
COMPUTE q51_bi_p_weight1=cbcl_q51_p1 * 0.489.
EXECUTE.
COMPUTE q56d_bi_p_weight1=cbcl_q56d_p1 * 0.342.
EXECUTE.
COMPUTE q56e_bi_p_weight1=cbcl_q56e_p1 * 0.311.
EXECUTE.
COMPUTE q111_bi_p_weight1=cbcl_q111_p1 * 0.653.
EXECUTE.
COMPUTE q42_bi_p_weight1=cbcl_q42_p1 * 0.512.
EXECUTE.
COMPUTE q75_bi_p_weight1=cbcl_q75_p1 * 0.418.
EXECUTE.
COMPUTE q65_bi_p_weight1=cbcl_q65_p1 * 0.602.
EXECUTE.
COMPUTE q102_bi_p_weight1=cbcl_q102_p1 * 0.591.
EXECUTE.

****Compute unstandardized factor scores for p from bifactor model

COMPUTE p_fs_bi1=SUM(comp1_bi_p_weight1,q16_bi_p_weight1,comp2_bi_p_weight1,q37_bi_p_weight1,
    q95_bi_p_weight1,q3_bi_p_weight1,comp3_bi_p_weight1,q68_bi_p_weight1,q26_bi_p_weight1,q90_bi_p_weight1,
    q94_bi_p_weight1,comp4_bi_p_weight1,q86_bi_p_weight1,q43_bi_p_weight1,q67_bi_p_weight1,q87_bi_p_weight1,
    q27_bi_p_weight1,comp5_bi_p_weight1,q89_bi_p_weight1,q19_bi_p_weight1,q39_bi_p_weight1,q34_bi_p_weight1,
    q88_bi_p_weight1,q7_bi_p_weight1,q109_bi_p_weight1,q50_bi_p_weight1,q112_bi_p_weight1,q32_bi_p_weight1,q52_bi_p_weight1,
    q45_bi_p_weight1,q31_bi_p_weight1,q35_bi_p_weight1,q30_bi_p_weight1,q29_bi_p_weight1,q12_bi_p_weight1,
    comp6_bi_p_weight1,q17_bi_p_weight1,q13_bi_p_weight1,q62_bi_p_weight1,
    q46_bi_p_weight1,q4_bi_p_weight1,q9_bi_p_weight1,q61_bi_p_weight1,q66_bi_p_weight1,q1_bi_p_weight1,
    q36_bi_p_weight1,q64_bi_p_weight1, q56c_bi_p_weight1,q56f_bi_p_weight1,q56g_bi_p_weight1,q56b_bi_p_weight1,
    q56a_bi_p_weight1,q51_bi_p_weight1,q56d_bi_p_weight1,q56e_bi_p_weight1,
    q111_bi_p_weight1,q42_bi_p_weight1,q75_bi_p_weight1,q65_bi_p_weight1,
    q102_bi_p_weight1).
EXECUTE.

****Wave 2*****

****EXT factor from higher-order model

COMPUTE comp1_weight2=comp1_2 * 0.441.
EXECUTE.
COMPUTE q16_weight2=cbcl_q16_p2 * 0.411.
EXECUTE.
COMPUTE comp2_weight2=comp2_2 * 0.456.
EXECUTE.
COMPUTE q37_weight2=cbcl_q37_p2 * 0.412.
EXECUTE.
COMPUTE q95_weight2=cbcl_q95_p2 * 0.444.
EXECUTE.
COMPUTE q3_weight2=cbcl_q03_p2 * 0.438.
EXECUTE.
COMPUTE comp3_weight2=comp3_2 * 0.438.
EXECUTE.
COMPUTE q68_weight2=cbcl_q68_p2 * 0.422.
EXECUTE.
COMPUTE q26_weight2=cbcl_q26_p2 * 0.401.
EXECUTE.
COMPUTE q90_weight2=cbcl_q90_p2 * 0.361.
EXECUTE.
COMPUTE q94_weight2=cbcl_q94_p2 * 0.379.
EXECUTE.
COMPUTE comp4_weight2=comp4_2 * 0.391.
EXECUTE.
COMPUTE q86_weight2=cbcl_q86_p2 * 0.440.
EXECUTE.
COMPUTE q43_weight2=cbcl_q43_p2 * 0.396.
EXECUTE.
COMPUTE q67_weight2=cbcl_q67_p2 * 0.410.
EXECUTE.
COMPUTE q87_weight2=cbcl_q87_p2 * 0.453.
EXECUTE.
COMPUTE q27_weight2=cbcl_q27_p2 * 0.395.
EXECUTE.
COMPUTE comp5_weight2=comp5_2 * 0.429.
EXECUTE.
COMPUTE q89_weight2=cbcl_q89_p2 * 0.420.
EXECUTE.
COMPUTE q19_weight2=cbcl_q19_p2 * 0.419.
EXECUTE.
COMPUTE q39_weight2=cbcl_q39_p2 * 0.314.
EXECUTE.
COMPUTE q34_weight2=cbcl_q34_p2 * 0.410.
EXECUTE.
COMPUTE q88_weight2=cbcl_q88_p2 * 0.417.
EXECUTE.
COMPUTE q7_weight2=cbcl_q07_p2 * 0.304.
EXECUTE.
COMPUTE q109_weight2=cbcl_q109_p2 * 0.361.
EXECUTE.

****Compute unstandardized factor scores for EXT from higher-order model

COMPUTE EXT_fs_ho2=SUM(comp1_weight2,q16_weight2,comp2_weight2,q37_weight2,q95_weight2,comp3_weight2,
    q68_weight2,q26_weight2,q90_weight2,q94_weight2,comp4_weight2,q86_weight2,q43_weight2,q67_weight2,
    q87_weight2,q27_weight2,comp5_weight2,q89_weight2,q19_weight2,q39_weight2,q34_weight2,q88_weight2,
    q109_weight2,q3_weight2,q7_weight2).
EXECUTE.

****INT factor from higher-order model

COMPUTE q50_weight2=cbcl_q50_p2 * 0.469.
EXECUTE.
COMPUTE q112_weight2=cbcl_q112_p2 * 0.441.
EXECUTE.
COMPUTE q32_weight2=cbcl_q32_p2 * 0.299.
EXECUTE.
COMPUTE q52_weight2=cbcl_q52_p2 * 0.425.
EXECUTE.
COMPUTE q45_weight2=cbcl_q45_p2 * 0.492.
EXECUTE.
COMPUTE q31_weight2=cbcl_q31_p2 * 0.403.
EXECUTE.
COMPUTE q35_weight2=cbcl_q35_p2 * 0.478.
EXECUTE.
COMPUTE q30_weight2=cbcl_q30_p2 * 0.410.
EXECUTE.
COMPUTE q29_weight2=cbcl_q29_p2 * 0.335.
EXECUTE.
COMPUTE q12_weight2=cbcl_q12_p2 * 0.455.
EXECUTE.

****Compute unstandardized factor scores for INT from higher-order model

COMPUTE INT_fs_ho2=SUM(q50_weight2,q112_weight2,q32_weight2,q52_weight2,q45_weight2,q31_weight2,q35_weight2,
    q30_weight2,q29_weight2,q12_weight2).
EXECUTE.

****ND factor from higher-order model

COMPUTE comp6_weight2=comp6_2 * 0.351.
EXECUTE.
COMPUTE q17_weight2=cbcl_q17_p2 * 0.266.
EXECUTE.
COMPUTE q13_weight2=cbcl_q13_p2 * 0.310.
EXECUTE.
COMPUTE q62_weight2=cbcl_q62_p2 * 0.292.
EXECUTE.
COMPUTE q46_weight2=cbcl_q46_p2 * 0.280.
EXECUTE.
COMPUTE q4_weight2=cbcl_q04_p2 * 0.336.
EXECUTE.
COMPUTE q9_weight2=cbcl_q09_p2 * 0.335.
EXECUTE.
COMPUTE q61_weight2=cbcl_q61_p2 * 0.308.
EXECUTE.
COMPUTE q66_weight2=cbcl_q66_p2 * 0.325.
EXECUTE.
COMPUTE q1_weight2=cbcl_q01_p2 * 0.280.
EXECUTE.
COMPUTE q36_weight2=cbcl_q36_p2 * 0.239.
EXECUTE.
COMPUTE q64_weight2=cbcl_q64_p2 * 0.246.
EXECUTE.

****Compute unstandardized factor scores for ND from higher-order model

COMPUTE ND_fs_ho2=SUM(comp6_weight2,q17_weight2,q13_weight2,q62_weight2,q46_weight2,q4_weight2,q9_weight2,
    q61_weight2,q66_weight2,q1_weight2,q36_weight2,q64_weight2).
EXECUTE.

****SOM factor from higher-order model

COMPUTE q56c_weight2=cbcl_q56c_p2 * 0.694.
EXECUTE.
COMPUTE q56f_weight2=cbcl_q56f_p2 * 0.623.
EXECUTE.
COMPUTE q56g_weight2=cbcl_q56g_p2 * 0.501.
EXECUTE.
COMPUTE q56b_weight2=cbcl_q56b_p2 * 0.527.
EXECUTE.
COMPUTE q56a_weight2=cbcl_q56a_p2 * 0.532.
EXECUTE.
COMPUTE q51_weight2=cbcl_q51_p2 * 0.641.
EXECUTE.
COMPUTE q56d_weight2=cbcl_q56d_p2 * 0.439.
EXECUTE.
COMPUTE q56e_weight2=cbcl_q56e_p2 * 0.403.
EXECUTE.

****Compute unstandardized factor scores for SOM from higher-order model

COMPUTE SOM_fs_ho2=SUM(q56c_weight2,q56f_weight2,q56g_weight2,q56b_weight2,q56a_weight2,q51_weight2,
    q56d_weight2,q56e_weight2).
EXECUTE.

****DET factor from higher-order model

COMPUTE q111_weight2=cbcl_q111_p2 * 0.566.
EXECUTE.
COMPUTE q42_weight2=cbcl_q42_p2 * 0.446.
EXECUTE.
COMPUTE q75_weight2=cbcl_q75_p2 * 0.366.
EXECUTE.
COMPUTE q65_weight2=cbcl_q65_p2 * 0.509.
EXECUTE.
COMPUTE q102_weight2=cbcl_q102_p2 * 0.491.
EXECUTE.

****Compute unstandardized factor scores for DET from higher-order model

COMPUTE DET_fs_ho2=SUM(q111_weight2,q42_weight2,q75_weight2,q65_weight2,q102_weight2).
EXECUTE.

****p factor from higher-order model

COMPUTE EXT_fs_ho_weight2=EXT_fs_ho3 * 1.465.
EXECUTE.
COMPUTE INT_fs_ho_weight2=INT_fs_ho3 * 1.357.
EXECUTE.
COMPUTE ND_fs_ho_weight2=ND_fs_ho3 * 2.069.
EXECUTE.
COMPUTE SOM_fs_ho_weight2=SOM_fs_ho3 * 0.691.
EXECUTE.
COMPUTE DET_fs_ho_weight2=DET_fs_ho3 * 1.169.
EXECUTE.

****Compute unstandardized factor scores for p from higher-order model

COMPUTE p_fs_ho2=SUM(EXT_fs_ho_weight2,INT_fs_ho_weight2,ND_fs_ho_weight2,SOM_fs_ho_weight2,
    DET_fs_ho_weight2).
EXECUTE.


****Multiply CBCL items by unstandardized weight2 from correlated factors model****

****EXT factor from correlated factors model

COMPUTE comp1_corr_weight2=comp1_2 * 0.782.
EXECUTE.
COMPUTE q16_corr_weight2=cbcl_q16_p2 * 0.729.
EXECUTE.
COMPUTE comp2_corr_weight2=comp2_2 * 0.809.
EXECUTE.
COMPUTE q37_corr_weight2=cbcl_q37_p2 * 0.731.
EXECUTE.
COMPUTE q95_corr_weight2=cbcl_q95_p2 * 0.787.
EXECUTE.
COMPUTE q3_corr_weight2=cbcl_q03_p2 * 0.776.
EXECUTE.
COMPUTE comp3_corr_weight2=comp3_2 * 0.778.
EXECUTE.
COMPUTE q68_corr_weight2=cbcl_q68_p2 * 0.749.
EXECUTE.
COMPUTE q26_corr_weight2=cbcl_q26_p2 * 0.711.
EXECUTE.
COMPUTE q90_corr_weight2=cbcl_q90_p2 * 0.641.
EXECUTE.
COMPUTE q94_corr_weight2=cbcl_q94_p2 * 0.672.
EXECUTE.
COMPUTE comp4_corr_weight2=comp4_2 * 0.694.
EXECUTE.
COMPUTE q86_corr_weight2=cbcl_q86_p2 * 0.781.
EXECUTE.
COMPUTE q43_corr_weight2=cbcl_q43_p2 * 0.703.
EXECUTE.
COMPUTE q67_corr_weight2=cbcl_q67_p2 * 0.728.
EXECUTE.
COMPUTE q87_corr_weight2=cbcl_q87_p2 * 0.803.
EXECUTE.
COMPUTE q27_corr_weight2=cbcl_q27_p2 * 0.701.
EXECUTE.
COMPUTE comp5_corr_weight2=comp5_2 * 0.760.
EXECUTE.
COMPUTE q89_corr_weight2=cbcl_q89_p2 * 0.744.
EXECUTE.
COMPUTE q19_corr_weight2=cbcl_q19_p2 * 0.744.
EXECUTE.
COMPUTE q39_corr_weight2=cbcl_q39_p2 * 0.558.
EXECUTE.
COMPUTE q34_corr_weight2=cbcl_q34_p2 * 0.726.
EXECUTE.
COMPUTE q88_corr_weight2=cbcl_q88_p2 * 0.739.
EXECUTE.
COMPUTE q7_corr_weight2=cbcl_q07_p2 * 0.539.
EXECUTE.
COMPUTE q109_corr_weight2=cbcl_q109_p2 * 0.640.
EXECUTE.

****Compute unstandardized factor scores for EXT from correlated factors model

COMPUTE EXT_fs_corr2=SUM(comp1_corr_weight2,q16_corr_weight2,comp2_corr_weight2,q37_corr_weight2,
    q95_corr_weight2,q3_corr_weight2,comp3_corr_weight2,q68_corr_weight2,q26_corr_weight2,q90_corr_weight2,
    q94_corr_weight2,comp4_corr_weight2,q86_corr_weight2,q43_corr_weight2,q67_corr_weight2,q87_corr_weight2,
    q27_corr_weight2,comp5_corr_weight2,q89_corr_weight2,q19_corr_weight2,q39_corr_weight2,q34_corr_weight2,
    q88_corr_weight2,q7_corr_weight2,q109_corr_weight2).
EXECUTE.

****INT factor from correlated factors model

COMPUTE q50_corr_weight2=cbcl_q50_p2 * 0.791.
EXECUTE.
COMPUTE q112_corr_weight2=cbcl_q112_p2 * 0.745.
EXECUTE.
COMPUTE q32_corr_weight2=cbcl_q32_p2 * 0.505.
EXECUTE.
COMPUTE q52_corr_weight2=cbcl_q52_p2 * 0.717.
EXECUTE.
COMPUTE q45_corr_weight2=cbcl_q45_p2 * 0.828.
EXECUTE.
COMPUTE q31_corr_weight2=cbcl_q31_p2 * 0.678.
EXECUTE.
COMPUTE q35_corr_weight2=cbcl_q35_p2 * 0.804.
EXECUTE.
COMPUTE q30_corr_weight2=cbcl_q30_p2 * 0.692.
EXECUTE.
COMPUTE q29_corr_weight2=cbcl_q29_p2 * 0.564.
EXECUTE.
COMPUTE q12_corr_weight2=cbcl_q12_p2 * 0.765.
EXECUTE.

****Compute unstandardized factor scores for INT from correlated factors model

COMPUTE INT_fs_corr2=SUM(q50_corr_weight2,q112_corr_weight2,q32_corr_weight2,q52_corr_weight2,
    q45_corr_weight2,q31_corr_weight2,q35_corr_weight2,q30_corr_weight2,q29_corr_weight2,q12_corr_weight2).
EXECUTE.

****ND factor from correlated factors model

COMPUTE comp6_corr_weight2=comp6_2 * 0.808.
EXECUTE.
COMPUTE q17_corr_weight2=cbcl_q17_p2 * 0.612.
EXECUTE.
COMPUTE q13_corr_weight2=cbcl_q13_p2 * 0.712.
EXECUTE.
COMPUTE q62_corr_weight2=cbcl_q62_p2 * 0.671.
EXECUTE.
COMPUTE q46_corr_weight2=cbcl_q46_p2 * 0.644.
EXECUTE.
COMPUTE q4_corr_weight2=cbcl_q04_p2 * 0.772.
EXECUTE.
COMPUTE q9_corr_weight2=cbcl_q09_p2 * 0.769.
EXECUTE.
COMPUTE q61_corr_weight2=cbcl_q61_p2 * 0.708.
EXECUTE.
COMPUTE q66_corr_weight2=cbcl_q66_p2 * 0.748.
EXECUTE.
COMPUTE q1_corr_weight2=cbcl_q01_p2 * 0.644.
EXECUTE.
COMPUTE q36_corr_weight2=cbcl_q36_p2 * 0.548.
EXECUTE.
COMPUTE q64_corr_weight2=cbcl_q64_p2 * 0.565.
EXECUTE.

****Compute unstandardized factor scores for ND from correlated factors model

COMPUTE ND_fs_corr2=SUM(comp6_corr_weight2,q17_corr_weight2,q13_corr_weight2,q62_corr_weight2,
    q46_corr_weight2,q4_corr_weight2,q9_corr_weight2,q61_corr_weight2,q66_corr_weight2,q1_corr_weight2,
    q36_corr_weight2,q64_corr_weight2).
EXECUTE.

****SOM factor from correlated factors model

COMPUTE q56c_corr_weight2=cbcl_q56c_p2 * 0.844.
EXECUTE.
COMPUTE q56f_corr_weight2=cbcl_q56f_p2 * 0.759.
EXECUTE.
COMPUTE q56g_corr_weight2=cbcl_q56g_p2 * 0.607.
EXECUTE.
COMPUTE q56b_corr_weight2=cbcl_q56b_p2 * 0.641.
EXECUTE.
COMPUTE q56a_corr_weight2=cbcl_q56a_p2 * 0.646.
EXECUTE.
COMPUTE q51_corr_weight2=cbcl_q51_p2 * 0.779.
EXECUTE.
COMPUTE q56d_corr_weight2=cbcl_q56d_p2 * 0.532.
EXECUTE.
COMPUTE q56e_corr_weight2=cbcl_q56e_p2 * 0.489.
EXECUTE.

****Compute unstandardized factor scores for SOM from correlated factors model

COMPUTE SOM_fs_corr2=SUM(q56c_corr_weight2,q56f_corr_weight2,q56g_corr_weight2,q56b_corr_weight2,
    q56a_corr_weight2,q51_corr_weight2,q56d_corr_weight2,q56e_corr_weight2).
EXECUTE.

****DET factor from correlated factors model

COMPUTE q111_corr_weight2=cbcl_q111_p2 * 0.869.
EXECUTE.
COMPUTE q42_corr_weight2=cbcl_q42_p2 * 0.686.
EXECUTE.
COMPUTE q75_corr_weight2=cbcl_q75_p2 * 0.567.
EXECUTE.
COMPUTE q65_corr_weight2=cbcl_q65_p2 * 0.781.
EXECUTE.
COMPUTE q102_corr_weight2=cbcl_q102_p2 * 0.754.
EXECUTE.

****Compute unstandardized factor scores for DET from correlated factors model

COMPUTE DET_fs_corr2=SUM(q111_corr_weight2,q42_corr_weight2,q75_corr_weight2,q65_corr_weight2,
    q102_corr_weight2).
EXECUTE.

****Multiply CBCL items by unstandardized weight2 from the bifactor model****

****EXT factor from bifactor model

COMPUTE comp1_bi_ext_weight2=comp1_2 * 0.628.
EXECUTE.
COMPUTE q16_bi_ext_weight2=cbcl_q16_p2 * 0.628.
EXECUTE.
COMPUTE comp2_bi_ext_weight2=comp2_2 * 0.558.
EXECUTE.
COMPUTE q37_bi_ext_weight2=cbcl_q37_p2 * 0.577.
EXECUTE.
COMPUTE q95_bi_ext_weight2=cbcl_q95_p2 * 0.473.
EXECUTE.
COMPUTE q3_bi_ext_weight2=cbcl_q03_p2 * 0.489.
EXECUTE.
COMPUTE comp3_bi_ext_weight2=comp3_2 * 0.484.
EXECUTE.
COMPUTE q68_bi_ext_weight2=cbcl_q68_p2 * 0.449.
EXECUTE.
COMPUTE q26_bi_ext_weight2=cbcl_q26_p2 * 0.502.
EXECUTE.
COMPUTE q90_bi_ext_weight2=cbcl_q90_p2 * 0.481.
EXECUTE.
COMPUTE q94_bi_ext_weight2=cbcl_q94_p2 * 0.493.
EXECUTE.
COMPUTE comp4_bi_ext_weight2=comp4_2 * 0.548.
EXECUTE.
COMPUTE q86_bi_ext_weight2=cbcl_q86_p2 * 0.364.
EXECUTE.
COMPUTE q43_bi_ext_weight2=cbcl_q43_p2 * 0.514.
EXECUTE.
COMPUTE q67_bi_ext_weight2=cbcl_q67_p2 * 0.325.
EXECUTE.
COMPUTE q87_bi_ext_weight2=cbcl_q87_p2 * 0.253.
EXECUTE.
COMPUTE q27_bi_ext_weight2=cbcl_q27_p2 * 0.312.
EXECUTE.
COMPUTE comp5_bi_ext_weight2=comp5_2 * 0.257.
EXECUTE.
COMPUTE q89_bi_ext_weight2=cbcl_q89_p2 * 0.267.
EXECUTE.
COMPUTE q19_bi_ext_weight2=cbcl_q19_p2 * 0.273.
EXECUTE.
COMPUTE q39_bi_ext_weight2=cbcl_q39_p2 * 0.405.
EXECUTE.
COMPUTE q34_bi_ext_weight2=cbcl_q34_p2 * 0.150.
EXECUTE.
COMPUTE q88_bi_ext_weight2=cbcl_q88_p2 * 0.156.
EXECUTE.
COMPUTE q7_bi_ext_weight2=cbcl_q07_p2 * 0.370.
EXECUTE.
COMPUTE q109_bi_ext_weight2=cbcl_q109_p2 * 0.189.
EXECUTE.

****Compute unstandardized factor scores for EXT from  bifactor model

COMPUTE EXT_fs_bi2=SUM(comp1_bi_ext_weight2,q16_bi_ext_weight2,comp2_bi_ext_weight2,q37_bi_ext_weight2,
    q95_bi_ext_weight2,q3_bi_ext_weight2,comp3_bi_ext_weight2,q68_bi_ext_weight2,q26_bi_ext_weight2,q90_bi_ext_weight2,
    q94_bi_ext_weight2,comp4_bi_ext_weight2,q86_bi_ext_weight2,q43_bi_ext_weight2,q67_bi_ext_weight2,q87_bi_ext_weight2,
    q27_bi_ext_weight2,comp5_bi_ext_weight2,q89_bi_ext_weight2,q19_bi_ext_weight2,q39_bi_ext_weight2,q34_bi_ext_weight2,
    q88_bi_ext_weight2,q7_bi_ext_weight2,q109_bi_ext_weight2).
EXECUTE.

****INT factor from bifactor model

COMPUTE q50_bi_int_weight2=cbcl_q50_p2 * 0.599.
EXECUTE.
COMPUTE q112_bi_int_weight2=cbcl_q112_p2 * 0.555.
EXECUTE.
COMPUTE q32_bi_int_weight2=cbcl_q32_p2 * 0.538.
EXECUTE.
COMPUTE q52_bi_int_weight2=cbcl_q52_p2 * 0.555.
EXECUTE.
COMPUTE q45_bi_int_weight2=cbcl_q45_p2 * 0.395.
EXECUTE.
COMPUTE q31_bi_int_weight2=cbcl_q31_p2 * 0.476.
EXECUTE.
COMPUTE q35_bi_int_weight2=cbcl_q35_p2 * 0.301.
EXECUTE.
COMPUTE q30_bi_int_weight2=cbcl_q30_p2 * 0.298.
EXECUTE.
COMPUTE q29_bi_int_weight2=cbcl_q29_p2 * 0.333.
EXECUTE.
COMPUTE q12_bi_int_weight2=cbcl_q12_p2 * 0.129.
EXECUTE.

****Compute unstandardized factor scores for INT from bifactor model

COMPUTE INT_fs_bi2=SUM(q50_bi_int_weight2,q112_bi_int_weight2,q32_bi_int_weight2,q52_bi_int_weight2,
    q45_bi_int_weight2,q31_bi_int_weight2,q35_bi_int_weight2,q30_bi_int_weight2,q29_bi_int_weight2,q12_bi_int_weight2).
EXECUTE.

****ND factor from bifactor model

COMPUTE comp6_bi_nd_weight2=comp6_2 * 0.489.
EXECUTE.
COMPUTE q17_bi_nd_weight2=cbcl_q17_p2 * 0.413.
EXECUTE.
COMPUTE q13_bi_nd_weight2=cbcl_q13_p2 * 0.383.
EXECUTE.
COMPUTE q62_bi_nd_weight2=cbcl_q62_p2 * 0.417.
EXECUTE.
COMPUTE q46_bi_nd_weight2=cbcl_q46_p2 * 0.144.
EXECUTE.
COMPUTE q4_bi_nd_weight2=cbcl_q04_p2 * 0.365.
EXECUTE.
COMPUTE q9_bi_nd_weight2=cbcl_q09_p2 * 0.129.
EXECUTE.
COMPUTE q61_bi_nd_weight2=cbcl_q61_p2 * 0.357.
EXECUTE.
COMPUTE q66_bi_nd_weight2=cbcl_q66_p2 * 0.137.
EXECUTE.
COMPUTE q1_bi_nd_weight2=cbcl_q01_p2 * 0.251.
EXECUTE.
COMPUTE q36_bi_nd_weight2=cbcl_q36_p2 * 0.293.
EXECUTE.
COMPUTE q64_bi_nd_weight2=cbcl_q64_p2 * 0.149.
EXECUTE.

****Compute unstandardized factor scores for ND from bifactor model

COMPUTE ND_fs_bi2=SUM(comp6_bi_nd_weight2,q17_bi_nd_weight2,q13_bi_nd_weight2,q62_bi_nd_weight2,
    q46_bi_nd_weight2,q4_bi_nd_weight2,q9_bi_nd_weight2,q61_bi_nd_weight2,q66_bi_nd_weight2,q1_bi_nd_weight2,
    q36_bi_nd_weight2,q64_bi_nd_weight2).
EXECUTE.

****SOM factor from bifactor model

COMPUTE q56c_bi_som_weight2=cbcl_q56c_p2 * 0.776.
EXECUTE.
COMPUTE q56f_bi_som_weight2=cbcl_q56f_p2 * 0.713.
EXECUTE.
COMPUTE q56g_bi_som_weight2=cbcl_q56g_p2 * 0.614.
EXECUTE.
COMPUTE q56b_bi_som_weight2=cbcl_q56b_p2 * 0.533.
EXECUTE.
COMPUTE q56a_bi_som_weight2=cbcl_q56a_p2 * 0.438.
EXECUTE.
COMPUTE q51_bi_som_weight2=cbcl_q51_p2 * 0.443.
EXECUTE.
COMPUTE q56d_bi_som_weight2=cbcl_q56d_p2 * 0.274.
EXECUTE.
COMPUTE q56e_bi_som_weight2=cbcl_q56e_p2 * 0.252.
EXECUTE.

****Compute unstandardized factor scores for SOM from bifactor model

COMPUTE SOM_fs_bi2=SUM(q56c_bi_som_weight2,q56f_bi_som_weight2,q56g_bi_som_weight2,q56b_bi_som_weight2,
    q56a_bi_som_weight2,q51_bi_som_weight2,q56d_bi_som_weight2,q56e_bi_som_weight2).
EXECUTE.

****DET factor from bifactor model

COMPUTE q111_bi_det_weight2=cbcl_q111_p2 * 0.630.
EXECUTE.
COMPUTE q42_bi_det_weight2=cbcl_q42_p2 * 0.531.
EXECUTE.
COMPUTE q75_bi_det_weight2=cbcl_q75_p2 * 0.474.
EXECUTE.
COMPUTE q65_bi_det_weight2=cbcl_q65_p2 * 0.398.
EXECUTE.
COMPUTE q102_bi_det_weight2=cbcl_q102_p2 * 0.299.
EXECUTE.

****Compute unstandardized factor scores for DET from bifactor model

COMPUTE DET_fs_bi2=SUM(q111_bi_det_weight2,q42_bi_det_weight2,q75_bi_det_weight2,q65_bi_det_weight2,
    q102_bi_det_weight2).
EXECUTE.

****p factor from bifactor model

COMPUTE comp1_bi_p_weight2=comp1_2 * 0.558.
EXECUTE.
COMPUTE q16_bi_p_weight2=cbcl_q16_p2 * 0.503.
EXECUTE.
COMPUTE comp2_bi_p_weight2=comp2_2 * 0.626.
EXECUTE.
COMPUTE q37_bi_p_weight2=cbcl_q37_p2 * 0.527.
EXECUTE.
COMPUTE q95_bi_p_weight2=cbcl_q95_p2 * 0.640.
EXECUTE.
COMPUTE q3_bi_p_weight2=cbcl_q03_p2 * 0.622.
EXECUTE.
COMPUTE comp3_bi_p_weight2=comp3_2 * 0.623.
EXECUTE.
COMPUTE q68_bi_p_weight2=cbcl_q68_p2 * 0.609.
EXECUTE.
COMPUTE q26_bi_p_weight2=cbcl_q26_p2 * 0.543.
EXECUTE.
COMPUTE q90_bi_p_weight2=cbcl_q90_p2 * 0.480.
EXECUTE.
COMPUTE q94_bi_p_weight2=cbcl_q94_p2 * 0.507.
EXECUTE.
COMPUTE comp4_bi_p_weight2=comp4_2 * 0.503.
EXECUTE.
COMPUTE q86_bi_p_weight2=cbcl_q86_p2 * 0.679.
EXECUTE.
COMPUTE q43_bi_p_weight2=cbcl_q43_p2 * 0.532.
EXECUTE.
COMPUTE q67_bi_p_weight2=cbcl_q67_p2 * 0.638.
EXECUTE.
COMPUTE q87_bi_p_weight2=cbcl_q87_p2 * 0.745.
EXECUTE.
COMPUTE q27_bi_p_weight2=cbcl_q27_p2 * 0.614.
EXECUTE.
COMPUTE comp5_bi_p_weight2=comp5_2 * 0.699.
EXECUTE.
COMPUTE q89_bi_p_weight2=cbcl_q89_p2 * 0.679.
EXECUTE.
COMPUTE q19_bi_p_weight2=cbcl_q19_p2 * 0.674.
EXECUTE.
COMPUTE q39_bi_p_weight2=cbcl_q39_p2 * 0.424.
EXECUTE.
COMPUTE q34_bi_p_weight2=cbcl_q34_p2 * 0.701.
EXECUTE.
COMPUTE q88_bi_p_weight2=cbcl_q88_p2 * 0.713.
EXECUTE.
COMPUTE q7_bi_p_weight2=cbcl_q07_p2 * 0.420.
EXECUTE.
COMPUTE q109_bi_p_weight2=cbcl_q109_p2 * 0.596.
EXECUTE.
COMPUTE q50_bi_p_weight2=cbcl_q50_p2 * 0.608.
EXECUTE.
COMPUTE q112_bi_p_weight2=cbcl_q112_p2 * 0.573.
EXECUTE.
COMPUTE q32_bi_p_weight2=cbcl_q32_p2 * 0.359.
EXECUTE.
COMPUTE q52_bi_p_weight2=cbcl_q52_p2 * 0.547.
EXECUTE.
COMPUTE q45_bi_p_weight2=cbcl_q45_p2 * 0.683.
EXECUTE.
COMPUTE q31_bi_p_weight2=cbcl_q31_p2 * 0.531.
EXECUTE.
COMPUTE q35_bi_p_weight2=cbcl_q35_p2 * 0.681.
EXECUTE.
COMPUTE q30_bi_p_weight2=cbcl_q30_p2 * 0.578.
EXECUTE.
COMPUTE q29_bi_p_weight2=cbcl_q29_p2 * 0.454.
EXECUTE.
COMPUTE q12_bi_p_weight2=cbcl_q12_p2 * 0.673.
EXECUTE.
COMPUTE comp6_bi_p_weight2=comp6_2 * 0.705.
EXECUTE.
COMPUTE q17_bi_p_weight2=cbcl_q17_p2 * 0.524.
EXECUTE.
COMPUTE q13_bi_p_weight2=cbcl_q13_p2 * 0.628.
EXECUTE.
COMPUTE q62_bi_p_weight2=cbcl_q62_p2 * 0.581.
EXECUTE.
COMPUTE q46_bi_p_weight2=cbcl_q46_p2 * 0.599.
EXECUTE.
COMPUTE q4_bi_p_weight2=cbcl_q04_p2 * 0.685.
EXECUTE.
COMPUTE q9_bi_p_weight2=cbcl_q09_p2 * 0.721.
EXECUTE.
COMPUTE q61_bi_p_weight2=cbcl_q61_p2 * 0.624.
EXECUTE.
COMPUTE q66_bi_p_weight2=cbcl_q66_p2 * 0.700.
EXECUTE.
COMPUTE q1_bi_p_weight2=cbcl_q01_p2 * 0.582.
EXECUTE.
COMPUTE q36_bi_p_weight2=cbcl_q36_p2 * 0.483.
EXECUTE.
COMPUTE q64_bi_p_weight2=cbcl_q64_p2 * 0.522.
EXECUTE.
COMPUTE q56c_bi_p_weight2=cbcl_q56c_p2 * 0.450.
EXECUTE.
COMPUTE q56f_bi_p_weight2=cbcl_q56f_p2 * 0.404.
EXECUTE.
COMPUTE q56g_bi_p_weight2=cbcl_q56g_p2 * 0.310.
EXECUTE.
COMPUTE q56b_bi_p_weight2=cbcl_q56b_p2 * 0.363.
EXECUTE.
COMPUTE q56a_bi_p_weight2=cbcl_q56a_p2 * 0.390.
EXECUTE.
COMPUTE q51_bi_p_weight2=cbcl_q51_p2 * 0.489.
EXECUTE.
COMPUTE q56d_bi_p_weight2=cbcl_q56d_p2 * 0.342.
EXECUTE.
COMPUTE q56e_bi_p_weight2=cbcl_q56e_p2 * 0.311.
EXECUTE.
COMPUTE q111_bi_p_weight2=cbcl_q111_p2 * 0.653.
EXECUTE.
COMPUTE q42_bi_p_weight2=cbcl_q42_p2 * 0.512.
EXECUTE.
COMPUTE q75_bi_p_weight2=cbcl_q75_p2 * 0.418.
EXECUTE.
COMPUTE q65_bi_p_weight2=cbcl_q65_p2 * 0.602.
EXECUTE.
COMPUTE q102_bi_p_weight2=cbcl_q102_p2 * 0.591.
EXECUTE.

****Compute unstandardized factor scores for p from bifactor model

COMPUTE p_fs_bi2=SUM(comp1_bi_p_weight2,q16_bi_p_weight2,comp2_bi_p_weight2,q37_bi_p_weight2,
    q95_bi_p_weight2,q3_bi_p_weight2,comp3_bi_p_weight2,q68_bi_p_weight2,q26_bi_p_weight2,q90_bi_p_weight2,
    q94_bi_p_weight2,comp4_bi_p_weight2,q86_bi_p_weight2,q43_bi_p_weight2,q67_bi_p_weight2,q87_bi_p_weight2,
    q27_bi_p_weight2,comp5_bi_p_weight2,q89_bi_p_weight2,q19_bi_p_weight2,q39_bi_p_weight2,q34_bi_p_weight2,
    q88_bi_p_weight2,q7_bi_p_weight2,q109_bi_p_weight2,q50_bi_p_weight2,q112_bi_p_weight2,q32_bi_p_weight2,q52_bi_p_weight2,
    q45_bi_p_weight2,q31_bi_p_weight2,q35_bi_p_weight2,q30_bi_p_weight2,q29_bi_p_weight2,q12_bi_p_weight2,
    comp6_bi_p_weight2,q17_bi_p_weight2,q13_bi_p_weight2,q62_bi_p_weight2,
    q46_bi_p_weight2,q4_bi_p_weight2,q9_bi_p_weight2,q61_bi_p_weight2,q66_bi_p_weight2,q1_bi_p_weight2,
    q36_bi_p_weight2,q64_bi_p_weight2, q56c_bi_p_weight2,q56f_bi_p_weight2,q56g_bi_p_weight2,q56b_bi_p_weight2,
    q56a_bi_p_weight2,q51_bi_p_weight2,q56d_bi_p_weight2,q56e_bi_p_weight2,
    q111_bi_p_weight2,q42_bi_p_weight2,q75_bi_p_weight2,q65_bi_p_weight2,
    q102_bi_p_weight2).
EXECUTE.

****Wave 3*****

****EXT factor from higher-order model

COMPUTE comp1_weight3=comp1_3 * 0.441.
EXECUTE.
COMPUTE q16_weight3=cbcl_q16_p3 * 0.411.
EXECUTE.
COMPUTE comp2_weight3=comp2_3 * 0.456.
EXECUTE.
COMPUTE q37_weight3=cbcl_q37_p3 * 0.412.
EXECUTE.
COMPUTE q95_weight3=cbcl_q95_p3 * 0.444.
EXECUTE.
COMPUTE q3_weight3=cbcl_q03_p3 * 0.438.
EXECUTE.
COMPUTE comp3_weight3=comp3_3 * 0.438.
EXECUTE.
COMPUTE q68_weight3=cbcl_q68_p3 * 0.422.
EXECUTE.
COMPUTE q26_weight3=cbcl_q26_p3 * 0.401.
EXECUTE.
COMPUTE q90_weight3=cbcl_q90_p3 * 0.361.
EXECUTE.
COMPUTE q94_weight3=cbcl_q94_p3 * 0.379.
EXECUTE.
COMPUTE comp4_weight3=comp4_3 * 0.391.
EXECUTE.
COMPUTE q86_weight3=cbcl_q86_p3 * 0.440.
EXECUTE.
COMPUTE q43_weight3=cbcl_q43_p3 * 0.396.
EXECUTE.
COMPUTE q67_weight3=cbcl_q67_p3 * 0.410.
EXECUTE.
COMPUTE q87_weight3=cbcl_q87_p3 * 0.453.
EXECUTE.
COMPUTE q27_weight3=cbcl_q27_p3 * 0.395.
EXECUTE.
COMPUTE comp5_weight3=comp5_3 * 0.429.
EXECUTE.
COMPUTE q89_weight3=cbcl_q89_p3 * 0.420.
EXECUTE.
COMPUTE q19_weight3=cbcl_q19_p3 * 0.419.
EXECUTE.
COMPUTE q39_weight3=cbcl_q39_p3 * 0.314.
EXECUTE.
COMPUTE q34_weight3=cbcl_q34_p3 * 0.410.
EXECUTE.
COMPUTE q88_weight3=cbcl_q88_p3 * 0.417.
EXECUTE.
COMPUTE q7_weight3=cbcl_q07_p3 * 0.304.
EXECUTE.
COMPUTE q109_weight3=cbcl_q109_p3 * 0.361.
EXECUTE.

****Compute unstandardized factor scores for EXT from higher-order model

COMPUTE EXT_fs_ho3=SUM(comp1_weight3,q16_weight3,comp2_weight3,q37_weight3,q95_weight3,comp3_weight3,
    q68_weight3,q26_weight3,q90_weight3,q94_weight3,comp4_weight3,q86_weight3,q43_weight3,q67_weight3,
    q87_weight3,q27_weight3,comp5_weight3,q89_weight3,q19_weight3,q39_weight3,q34_weight3,q88_weight3,
    q109_weight3,q3_weight3,q7_weight3).
EXECUTE.

****INT factor from higher-order model

COMPUTE q50_weight3=cbcl_q50_p3 * 0.469.
EXECUTE.
COMPUTE q112_weight3=cbcl_q112_p3 * 0.441.
EXECUTE.
COMPUTE q32_weight3=cbcl_q32_p3 * 0.299.
EXECUTE.
COMPUTE q52_weight3=cbcl_q52_p3 * 0.425.
EXECUTE.
COMPUTE q45_weight3=cbcl_q45_p3 * 0.492.
EXECUTE.
COMPUTE q31_weight3=cbcl_q31_p3 * 0.403.
EXECUTE.
COMPUTE q35_weight3=cbcl_q35_p3 * 0.478.
EXECUTE.
COMPUTE q30_weight3=cbcl_q30_p3 * 0.410.
EXECUTE.
COMPUTE q29_weight3=cbcl_q29_p3 * 0.335.
EXECUTE.
COMPUTE q12_weight3=cbcl_q12_p3 * 0.455.
EXECUTE.

****Compute unstandardized factor scores for INT from higher-order model

COMPUTE INT_fs_ho3=SUM(q50_weight3,q112_weight3,q32_weight3,q52_weight3,q45_weight3,q31_weight3,q35_weight3,
    q30_weight3,q29_weight3,q12_weight3).
EXECUTE.

****ND factor from higher-order model

COMPUTE comp6_weight3=comp6_3 * 0.351.
EXECUTE.
COMPUTE q17_weight3=cbcl_q17_p3 * 0.266.
EXECUTE.
COMPUTE q13_weight3=cbcl_q13_p3 * 0.310.
EXECUTE.
COMPUTE q62_weight3=cbcl_q62_p3 * 0.292.
EXECUTE.
COMPUTE q46_weight3=cbcl_q46_p3 * 0.280.
EXECUTE.
COMPUTE q4_weight3=cbcl_q04_p3 * 0.336.
EXECUTE.
COMPUTE q9_weight3=cbcl_q09_p3 * 0.335.
EXECUTE.
COMPUTE q61_weight3=cbcl_q61_p3 * 0.308.
EXECUTE.
COMPUTE q66_weight3=cbcl_q66_p3 * 0.325.
EXECUTE.
COMPUTE q1_weight3=cbcl_q01_p3 * 0.280.
EXECUTE.
COMPUTE q36_weight3=cbcl_q36_p3 * 0.239.
EXECUTE.
COMPUTE q64_weight3=cbcl_q64_p3 * 0.246.
EXECUTE.

****Compute unstandardized factor scores for ND from higher-order model

COMPUTE ND_fs_ho3=SUM(comp6_weight3,q17_weight3,q13_weight3,q62_weight3,q46_weight3,q4_weight3,q9_weight3,
    q61_weight3,q66_weight3,q1_weight3,q36_weight3,q64_weight3).
EXECUTE.

****SOM factor from higher-order model

COMPUTE q56c_weight3=cbcl_q56c_p3 * 0.694.
EXECUTE.
COMPUTE q56f_weight3=cbcl_q56f_p3 * 0.623.
EXECUTE.
COMPUTE q56g_weight3=cbcl_q56g_p3 * 0.501.
EXECUTE.
COMPUTE q56b_weight3=cbcl_q56b_p3 * 0.527.
EXECUTE.
COMPUTE q56a_weight3=cbcl_q56a_p3 * 0.532.
EXECUTE.
COMPUTE q51_weight3=cbcl_q51_p3 * 0.641.
EXECUTE.
COMPUTE q56d_weight3=cbcl_q56d_p3 * 0.439.
EXECUTE.
COMPUTE q56e_weight3=cbcl_q56e_p3 * 0.403.
EXECUTE.

****Compute unstandardized factor scores for SOM from higher-order model

COMPUTE SOM_fs_ho3=SUM(q56c_weight3,q56f_weight3,q56g_weight3,q56b_weight3,q56a_weight3,q51_weight3,
    q56d_weight3,q56e_weight3).
EXECUTE.

****DET factor from higher-order model

COMPUTE q111_weight3=cbcl_q111_p3 * 0.566.
EXECUTE.
COMPUTE q42_weight3=cbcl_q42_p3 * 0.446.
EXECUTE.
COMPUTE q75_weight3=cbcl_q75_p3 * 0.366.
EXECUTE.
COMPUTE q65_weight3=cbcl_q65_p3 * 0.509.
EXECUTE.
COMPUTE q102_weight3=cbcl_q102_p3 * 0.491.
EXECUTE.

****Compute unstandardized factor scores for DET from higher-order model

COMPUTE DET_fs_ho3=SUM(q111_weight3,q42_weight3,q75_weight3,q65_weight3,q102_weight3).
EXECUTE.

****p factor from higher-order model

COMPUTE EXT_fs_ho_weight3=EXT_fs_ho3 * 1.465.
EXECUTE.
COMPUTE INT_fs_ho_weight3=INT_fs_ho3 * 1.357.
EXECUTE.
COMPUTE ND_fs_ho_weight3=ND_fs_ho3 * 2.069.
EXECUTE.
COMPUTE SOM_fs_ho_weight3=SOM_fs_ho3 * 0.691.
EXECUTE.
COMPUTE DET_fs_ho_weight3=DET_fs_ho3 * 1.169.
EXECUTE.

****Compute unstandardized factor scores for p from higher-order model

COMPUTE p_fs_ho3=SUM(EXT_fs_ho_weight3,INT_fs_ho_weight3,ND_fs_ho_weight3,SOM_fs_ho_weight3,
    DET_fs_ho_weight3).
EXECUTE.


****Multiply CBCL items by unstandardized weight3 from correlated factors model****

****EXT factor from correlated factors model

COMPUTE comp1_corr_weight3=comp1_3 * 0.782.
EXECUTE.
COMPUTE q16_corr_weight3=cbcl_q16_p3 * 0.729.
EXECUTE.
COMPUTE comp2_corr_weight3=comp2_3 * 0.809.
EXECUTE.
COMPUTE q37_corr_weight3=cbcl_q37_p3 * 0.731.
EXECUTE.
COMPUTE q95_corr_weight3=cbcl_q95_p3 * 0.787.
EXECUTE.
COMPUTE q3_corr_weight3=cbcl_q03_p3 * 0.776.
EXECUTE.
COMPUTE comp3_corr_weight3=comp3_3 * 0.778.
EXECUTE.
COMPUTE q68_corr_weight3=cbcl_q68_p3 * 0.749.
EXECUTE.
COMPUTE q26_corr_weight3=cbcl_q26_p3 * 0.711.
EXECUTE.
COMPUTE q90_corr_weight3=cbcl_q90_p3 * 0.641.
EXECUTE.
COMPUTE q94_corr_weight3=cbcl_q94_p3 * 0.672.
EXECUTE.
COMPUTE comp4_corr_weight3=comp4_3 * 0.694.
EXECUTE.
COMPUTE q86_corr_weight3=cbcl_q86_p3 * 0.781.
EXECUTE.
COMPUTE q43_corr_weight3=cbcl_q43_p3 * 0.703.
EXECUTE.
COMPUTE q67_corr_weight3=cbcl_q67_p3 * 0.728.
EXECUTE.
COMPUTE q87_corr_weight3=cbcl_q87_p3 * 0.803.
EXECUTE.
COMPUTE q27_corr_weight3=cbcl_q27_p3 * 0.701.
EXECUTE.
COMPUTE comp5_corr_weight3=comp5_3 * 0.760.
EXECUTE.
COMPUTE q89_corr_weight3=cbcl_q89_p3 * 0.744.
EXECUTE.
COMPUTE q19_corr_weight3=cbcl_q19_p3 * 0.744.
EXECUTE.
COMPUTE q39_corr_weight3=cbcl_q39_p3 * 0.558.
EXECUTE.
COMPUTE q34_corr_weight3=cbcl_q34_p3 * 0.726.
EXECUTE.
COMPUTE q88_corr_weight3=cbcl_q88_p3 * 0.739.
EXECUTE.
COMPUTE q7_corr_weight3=cbcl_q07_p3 * 0.539.
EXECUTE.
COMPUTE q109_corr_weight3=cbcl_q109_p3 * 0.640.
EXECUTE.

****Compute unstandardized factor scores for EXT from correlated factors model

COMPUTE EXT_fs_corr3=SUM(comp1_corr_weight3,q16_corr_weight3,comp2_corr_weight3,q37_corr_weight3,
    q95_corr_weight3,q3_corr_weight3,comp3_corr_weight3,q68_corr_weight3,q26_corr_weight3,q90_corr_weight3,
    q94_corr_weight3,comp4_corr_weight3,q86_corr_weight3,q43_corr_weight3,q67_corr_weight3,q87_corr_weight3,
    q27_corr_weight3,comp5_corr_weight3,q89_corr_weight3,q19_corr_weight3,q39_corr_weight3,q34_corr_weight3,
    q88_corr_weight3,q7_corr_weight3,q109_corr_weight3).
EXECUTE.

****INT factor from correlated factors model

COMPUTE q50_corr_weight3=cbcl_q50_p3 * 0.791.
EXECUTE.
COMPUTE q112_corr_weight3=cbcl_q112_p3 * 0.745.
EXECUTE.
COMPUTE q32_corr_weight3=cbcl_q32_p3 * 0.505.
EXECUTE.
COMPUTE q52_corr_weight3=cbcl_q52_p3 * 0.717.
EXECUTE.
COMPUTE q45_corr_weight3=cbcl_q45_p3 * 0.828.
EXECUTE.
COMPUTE q31_corr_weight3=cbcl_q31_p3 * 0.678.
EXECUTE.
COMPUTE q35_corr_weight3=cbcl_q35_p3 * 0.804.
EXECUTE.
COMPUTE q30_corr_weight3=cbcl_q30_p3 * 0.692.
EXECUTE.
COMPUTE q29_corr_weight3=cbcl_q29_p3 * 0.564.
EXECUTE.
COMPUTE q12_corr_weight3=cbcl_q12_p3 * 0.765.
EXECUTE.

****Compute unstandardized factor scores for INT from correlated factors model

COMPUTE INT_fs_corr3=SUM(q50_corr_weight3,q112_corr_weight3,q32_corr_weight3,q52_corr_weight3,
    q45_corr_weight3,q31_corr_weight3,q35_corr_weight3,q30_corr_weight3,q29_corr_weight3,q12_corr_weight3).
EXECUTE.

****ND factor from correlated factors model

COMPUTE comp6_corr_weight3=comp6_3 * 0.808.
EXECUTE.
COMPUTE q17_corr_weight3=cbcl_q17_p3 * 0.612.
EXECUTE.
COMPUTE q13_corr_weight3=cbcl_q13_p3 * 0.712.
EXECUTE.
COMPUTE q62_corr_weight3=cbcl_q62_p3 * 0.671.
EXECUTE.
COMPUTE q46_corr_weight3=cbcl_q46_p3 * 0.644.
EXECUTE.
COMPUTE q4_corr_weight3=cbcl_q04_p3 * 0.772.
EXECUTE.
COMPUTE q9_corr_weight3=cbcl_q09_p3 * 0.769.
EXECUTE.
COMPUTE q61_corr_weight3=cbcl_q61_p3 * 0.708.
EXECUTE.
COMPUTE q66_corr_weight3=cbcl_q66_p3 * 0.748.
EXECUTE.
COMPUTE q1_corr_weight3=cbcl_q01_p3 * 0.644.
EXECUTE.
COMPUTE q36_corr_weight3=cbcl_q36_p3 * 0.548.
EXECUTE.
COMPUTE q64_corr_weight3=cbcl_q64_p3 * 0.565.
EXECUTE.

****Compute unstandardized factor scores for ND from correlated factors model

COMPUTE ND_fs_corr3=SUM(comp6_corr_weight3,q17_corr_weight3,q13_corr_weight3,q62_corr_weight3,
    q46_corr_weight3,q4_corr_weight3,q9_corr_weight3,q61_corr_weight3,q66_corr_weight3,q1_corr_weight3,
    q36_corr_weight3,q64_corr_weight3).
EXECUTE.

****SOM factor from correlated factors model

COMPUTE q56c_corr_weight3=cbcl_q56c_p3 * 0.844.
EXECUTE.
COMPUTE q56f_corr_weight3=cbcl_q56f_p3 * 0.759.
EXECUTE.
COMPUTE q56g_corr_weight3=cbcl_q56g_p3 * 0.607.
EXECUTE.
COMPUTE q56b_corr_weight3=cbcl_q56b_p3 * 0.641.
EXECUTE.
COMPUTE q56a_corr_weight3=cbcl_q56a_p3 * 0.646.
EXECUTE.
COMPUTE q51_corr_weight3=cbcl_q51_p3 * 0.779.
EXECUTE.
COMPUTE q56d_corr_weight3=cbcl_q56d_p3 * 0.532.
EXECUTE.
COMPUTE q56e_corr_weight3=cbcl_q56e_p3 * 0.489.
EXECUTE.

****Compute unstandardized factor scores for SOM from correlated factors model

COMPUTE SOM_fs_corr3=SUM(q56c_corr_weight3,q56f_corr_weight3,q56g_corr_weight3,q56b_corr_weight3,
    q56a_corr_weight3,q51_corr_weight3,q56d_corr_weight3,q56e_corr_weight3).
EXECUTE.

****DET factor from correlated factors model

COMPUTE q111_corr_weight3=cbcl_q111_p3 * 0.869.
EXECUTE.
COMPUTE q42_corr_weight3=cbcl_q42_p3 * 0.686.
EXECUTE.
COMPUTE q75_corr_weight3=cbcl_q75_p3 * 0.567.
EXECUTE.
COMPUTE q65_corr_weight3=cbcl_q65_p3 * 0.781.
EXECUTE.
COMPUTE q102_corr_weight3=cbcl_q102_p3 * 0.754.
EXECUTE.

****Compute unstandardized factor scores for DET from correlated factors model

COMPUTE DET_fs_corr3=SUM(q111_corr_weight3,q42_corr_weight3,q75_corr_weight3,q65_corr_weight3,
    q102_corr_weight3).
EXECUTE.

****Multiply CBCL items by unstandardized weight3 from the bifactor model****

****EXT factor from bifactor model

COMPUTE comp1_bi_ext_weight3=comp1_3 * 0.628.
EXECUTE.
COMPUTE q16_bi_ext_weight3=cbcl_q16_p3 * 0.628.
EXECUTE.
COMPUTE comp2_bi_ext_weight3=comp2_3 * 0.558.
EXECUTE.
COMPUTE q37_bi_ext_weight3=cbcl_q37_p3 * 0.577.
EXECUTE.
COMPUTE q95_bi_ext_weight3=cbcl_q95_p3 * 0.473.
EXECUTE.
COMPUTE q3_bi_ext_weight3=cbcl_q03_p3 * 0.489.
EXECUTE.
COMPUTE comp3_bi_ext_weight3=comp3_3 * 0.484.
EXECUTE.
COMPUTE q68_bi_ext_weight3=cbcl_q68_p3 * 0.449.
EXECUTE.
COMPUTE q26_bi_ext_weight3=cbcl_q26_p3 * 0.502.
EXECUTE.
COMPUTE q90_bi_ext_weight3=cbcl_q90_p3 * 0.481.
EXECUTE.
COMPUTE q94_bi_ext_weight3=cbcl_q94_p3 * 0.493.
EXECUTE.
COMPUTE comp4_bi_ext_weight3=comp4_3 * 0.548.
EXECUTE.
COMPUTE q86_bi_ext_weight3=cbcl_q86_p3 * 0.364.
EXECUTE.
COMPUTE q43_bi_ext_weight3=cbcl_q43_p3 * 0.514.
EXECUTE.
COMPUTE q67_bi_ext_weight3=cbcl_q67_p3 * 0.325.
EXECUTE.
COMPUTE q87_bi_ext_weight3=cbcl_q87_p3 * 0.253.
EXECUTE.
COMPUTE q27_bi_ext_weight3=cbcl_q27_p3 * 0.312.
EXECUTE.
COMPUTE comp5_bi_ext_weight3=comp5_3 * 0.257.
EXECUTE.
COMPUTE q89_bi_ext_weight3=cbcl_q89_p3 * 0.267.
EXECUTE.
COMPUTE q19_bi_ext_weight3=cbcl_q19_p3 * 0.273.
EXECUTE.
COMPUTE q39_bi_ext_weight3=cbcl_q39_p3 * 0.405.
EXECUTE.
COMPUTE q34_bi_ext_weight3=cbcl_q34_p3 * 0.150.
EXECUTE.
COMPUTE q88_bi_ext_weight3=cbcl_q88_p3 * 0.156.
EXECUTE.
COMPUTE q7_bi_ext_weight3=cbcl_q07_p3 * 0.370.
EXECUTE.
COMPUTE q109_bi_ext_weight3=cbcl_q109_p3 * 0.189.
EXECUTE.

****Compute unstandardized factor scores for EXT from  bifactor model

COMPUTE EXT_fs_bi3=SUM(comp1_bi_ext_weight3,q16_bi_ext_weight3,comp2_bi_ext_weight3,q37_bi_ext_weight3,
    q95_bi_ext_weight3,q3_bi_ext_weight3,comp3_bi_ext_weight3,q68_bi_ext_weight3,q26_bi_ext_weight3,q90_bi_ext_weight3,
    q94_bi_ext_weight3,comp4_bi_ext_weight3,q86_bi_ext_weight3,q43_bi_ext_weight3,q67_bi_ext_weight3,q87_bi_ext_weight3,
    q27_bi_ext_weight3,comp5_bi_ext_weight3,q89_bi_ext_weight3,q19_bi_ext_weight3,q39_bi_ext_weight3,q34_bi_ext_weight3,
    q88_bi_ext_weight3,q7_bi_ext_weight3,q109_bi_ext_weight3).
EXECUTE.

****INT factor from bifactor model

COMPUTE q50_bi_int_weight3=cbcl_q50_p3 * 0.599.
EXECUTE.
COMPUTE q112_bi_int_weight3=cbcl_q112_p3 * 0.555.
EXECUTE.
COMPUTE q32_bi_int_weight3=cbcl_q32_p3 * 0.538.
EXECUTE.
COMPUTE q52_bi_int_weight3=cbcl_q52_p3 * 0.555.
EXECUTE.
COMPUTE q45_bi_int_weight3=cbcl_q45_p3 * 0.395.
EXECUTE.
COMPUTE q31_bi_int_weight3=cbcl_q31_p3 * 0.476.
EXECUTE.
COMPUTE q35_bi_int_weight3=cbcl_q35_p3 * 0.301.
EXECUTE.
COMPUTE q30_bi_int_weight3=cbcl_q30_p3 * 0.298.
EXECUTE.
COMPUTE q29_bi_int_weight3=cbcl_q29_p3 * 0.333.
EXECUTE.
COMPUTE q12_bi_int_weight3=cbcl_q12_p3 * 0.129.
EXECUTE.

****Compute unstandardized factor scores for INT from bifactor model

COMPUTE INT_fs_bi3=SUM(q50_bi_int_weight3,q112_bi_int_weight3,q32_bi_int_weight3,q52_bi_int_weight3,
    q45_bi_int_weight3,q31_bi_int_weight3,q35_bi_int_weight3,q30_bi_int_weight3,q29_bi_int_weight3,q12_bi_int_weight3).
EXECUTE.

****ND factor from bifactor model

COMPUTE comp6_bi_nd_weight3=comp6_3 * 0.489.
EXECUTE.
COMPUTE q17_bi_nd_weight3=cbcl_q17_p3 * 0.413.
EXECUTE.
COMPUTE q13_bi_nd_weight3=cbcl_q13_p3 * 0.383.
EXECUTE.
COMPUTE q62_bi_nd_weight3=cbcl_q62_p3 * 0.417.
EXECUTE.
COMPUTE q46_bi_nd_weight3=cbcl_q46_p3 * 0.144.
EXECUTE.
COMPUTE q4_bi_nd_weight3=cbcl_q04_p3 * 0.365.
EXECUTE.
COMPUTE q9_bi_nd_weight3=cbcl_q09_p3 * 0.129.
EXECUTE.
COMPUTE q61_bi_nd_weight3=cbcl_q61_p3 * 0.357.
EXECUTE.
COMPUTE q66_bi_nd_weight3=cbcl_q66_p3 * 0.137.
EXECUTE.
COMPUTE q1_bi_nd_weight3=cbcl_q01_p3 * 0.251.
EXECUTE.
COMPUTE q36_bi_nd_weight3=cbcl_q36_p3 * 0.293.
EXECUTE.
COMPUTE q64_bi_nd_weight3=cbcl_q64_p3 * 0.149.
EXECUTE.

****Compute unstandardized factor scores for ND from bifactor model

COMPUTE ND_fs_bi3=SUM(comp6_bi_nd_weight3,q17_bi_nd_weight3,q13_bi_nd_weight3,q62_bi_nd_weight3,
    q46_bi_nd_weight3,q4_bi_nd_weight3,q9_bi_nd_weight3,q61_bi_nd_weight3,q66_bi_nd_weight3,q1_bi_nd_weight3,
    q36_bi_nd_weight3,q64_bi_nd_weight3).
EXECUTE.

****SOM factor from bifactor model

COMPUTE q56c_bi_som_weight3=cbcl_q56c_p3 * 0.776.
EXECUTE.
COMPUTE q56f_bi_som_weight3=cbcl_q56f_p3 * 0.713.
EXECUTE.
COMPUTE q56g_bi_som_weight3=cbcl_q56g_p3 * 0.614.
EXECUTE.
COMPUTE q56b_bi_som_weight3=cbcl_q56b_p3 * 0.533.
EXECUTE.
COMPUTE q56a_bi_som_weight3=cbcl_q56a_p3 * 0.438.
EXECUTE.
COMPUTE q51_bi_som_weight3=cbcl_q51_p3 * 0.443.
EXECUTE.
COMPUTE q56d_bi_som_weight3=cbcl_q56d_p3 * 0.274.
EXECUTE.
COMPUTE q56e_bi_som_weight3=cbcl_q56e_p3 * 0.252.
EXECUTE.

****Compute unstandardized factor scores for SOM from bifactor model

COMPUTE SOM_fs_bi3=SUM(q56c_bi_som_weight3,q56f_bi_som_weight3,q56g_bi_som_weight3,q56b_bi_som_weight3,
    q56a_bi_som_weight3,q51_bi_som_weight3,q56d_bi_som_weight3,q56e_bi_som_weight3).
EXECUTE.

****DET factor from bifactor model

COMPUTE q111_bi_det_weight3=cbcl_q111_p3 * 0.630.
EXECUTE.
COMPUTE q42_bi_det_weight3=cbcl_q42_p3 * 0.531.
EXECUTE.
COMPUTE q75_bi_det_weight3=cbcl_q75_p3 * 0.474.
EXECUTE.
COMPUTE q65_bi_det_weight3=cbcl_q65_p3 * 0.398.
EXECUTE.
COMPUTE q102_bi_det_weight3=cbcl_q102_p3 * 0.299.
EXECUTE.

****Compute unstandardized factor scores for DET from bifactor model

COMPUTE DET_fs_bi3=SUM(q111_bi_det_weight3,q42_bi_det_weight3,q75_bi_det_weight3,q65_bi_det_weight3,
    q102_bi_det_weight3).
EXECUTE.

****p factor from bifactor model

COMPUTE comp1_bi_p_weight3=comp1_3 * 0.558.
EXECUTE.
COMPUTE q16_bi_p_weight3=cbcl_q16_p3 * 0.503.
EXECUTE.
COMPUTE comp2_bi_p_weight3=comp2_3 * 0.626.
EXECUTE.
COMPUTE q37_bi_p_weight3=cbcl_q37_p3 * 0.527.
EXECUTE.
COMPUTE q95_bi_p_weight3=cbcl_q95_p3 * 0.640.
EXECUTE.
COMPUTE q3_bi_p_weight3=cbcl_q03_p3 * 0.622.
EXECUTE.
COMPUTE comp3_bi_p_weight3=comp3_3 * 0.623.
EXECUTE.
COMPUTE q68_bi_p_weight3=cbcl_q68_p3 * 0.609.
EXECUTE.
COMPUTE q26_bi_p_weight3=cbcl_q26_p3 * 0.543.
EXECUTE.
COMPUTE q90_bi_p_weight3=cbcl_q90_p3 * 0.480.
EXECUTE.
COMPUTE q94_bi_p_weight3=cbcl_q94_p3 * 0.507.
EXECUTE.
COMPUTE comp4_bi_p_weight3=comp4_3 * 0.503.
EXECUTE.
COMPUTE q86_bi_p_weight3=cbcl_q86_p3 * 0.679.
EXECUTE.
COMPUTE q43_bi_p_weight3=cbcl_q43_p3 * 0.532.
EXECUTE.
COMPUTE q67_bi_p_weight3=cbcl_q67_p3 * 0.638.
EXECUTE.
COMPUTE q87_bi_p_weight3=cbcl_q87_p3 * 0.745.
EXECUTE.
COMPUTE q27_bi_p_weight3=cbcl_q27_p3 * 0.614.
EXECUTE.
COMPUTE comp5_bi_p_weight3=comp5_3 * 0.699.
EXECUTE.
COMPUTE q89_bi_p_weight3=cbcl_q89_p3 * 0.679.
EXECUTE.
COMPUTE q19_bi_p_weight3=cbcl_q19_p3 * 0.674.
EXECUTE.
COMPUTE q39_bi_p_weight3=cbcl_q39_p3 * 0.424.
EXECUTE.
COMPUTE q34_bi_p_weight3=cbcl_q34_p3 * 0.701.
EXECUTE.
COMPUTE q88_bi_p_weight3=cbcl_q88_p3 * 0.713.
EXECUTE.
COMPUTE q7_bi_p_weight3=cbcl_q07_p3 * 0.420.
EXECUTE.
COMPUTE q109_bi_p_weight3=cbcl_q109_p3 * 0.596.
EXECUTE.
COMPUTE q50_bi_p_weight3=cbcl_q50_p3 * 0.608.
EXECUTE.
COMPUTE q112_bi_p_weight3=cbcl_q112_p3 * 0.573.
EXECUTE.
COMPUTE q32_bi_p_weight3=cbcl_q32_p3 * 0.359.
EXECUTE.
COMPUTE q52_bi_p_weight3=cbcl_q52_p3 * 0.547.
EXECUTE.
COMPUTE q45_bi_p_weight3=cbcl_q45_p3 * 0.683.
EXECUTE.
COMPUTE q31_bi_p_weight3=cbcl_q31_p3 * 0.531.
EXECUTE.
COMPUTE q35_bi_p_weight3=cbcl_q35_p3 * 0.681.
EXECUTE.
COMPUTE q30_bi_p_weight3=cbcl_q30_p3 * 0.578.
EXECUTE.
COMPUTE q29_bi_p_weight3=cbcl_q29_p3 * 0.454.
EXECUTE.
COMPUTE q12_bi_p_weight3=cbcl_q12_p3 * 0.673.
EXECUTE.
COMPUTE comp6_bi_p_weight3=comp6_3 * 0.705.
EXECUTE.
COMPUTE q17_bi_p_weight3=cbcl_q17_p3 * 0.524.
EXECUTE.
COMPUTE q13_bi_p_weight3=cbcl_q13_p3 * 0.628.
EXECUTE.
COMPUTE q62_bi_p_weight3=cbcl_q62_p3 * 0.581.
EXECUTE.
COMPUTE q46_bi_p_weight3=cbcl_q46_p3 * 0.599.
EXECUTE.
COMPUTE q4_bi_p_weight3=cbcl_q04_p3 * 0.685.
EXECUTE.
COMPUTE q9_bi_p_weight3=cbcl_q09_p3 * 0.721.
EXECUTE.
COMPUTE q61_bi_p_weight3=cbcl_q61_p3 * 0.624.
EXECUTE.
COMPUTE q66_bi_p_weight3=cbcl_q66_p3 * 0.700.
EXECUTE.
COMPUTE q1_bi_p_weight3=cbcl_q01_p3 * 0.582.
EXECUTE.
COMPUTE q36_bi_p_weight3=cbcl_q36_p3 * 0.483.
EXECUTE.
COMPUTE q64_bi_p_weight3=cbcl_q64_p3 * 0.522.
EXECUTE.
COMPUTE q56c_bi_p_weight3=cbcl_q56c_p3 * 0.450.
EXECUTE.
COMPUTE q56f_bi_p_weight3=cbcl_q56f_p3 * 0.404.
EXECUTE.
COMPUTE q56g_bi_p_weight3=cbcl_q56g_p3 * 0.310.
EXECUTE.
COMPUTE q56b_bi_p_weight3=cbcl_q56b_p3 * 0.363.
EXECUTE.
COMPUTE q56a_bi_p_weight3=cbcl_q56a_p3 * 0.390.
EXECUTE.
COMPUTE q51_bi_p_weight3=cbcl_q51_p3 * 0.489.
EXECUTE.
COMPUTE q56d_bi_p_weight3=cbcl_q56d_p3 * 0.342.
EXECUTE.
COMPUTE q56e_bi_p_weight3=cbcl_q56e_p3 * 0.311.
EXECUTE.
COMPUTE q111_bi_p_weight3=cbcl_q111_p3 * 0.653.
EXECUTE.
COMPUTE q42_bi_p_weight3=cbcl_q42_p3 * 0.512.
EXECUTE.
COMPUTE q75_bi_p_weight3=cbcl_q75_p3 * 0.418.
EXECUTE.
COMPUTE q65_bi_p_weight3=cbcl_q65_p3 * 0.602.
EXECUTE.
COMPUTE q102_bi_p_weight3=cbcl_q102_p3 * 0.591.
EXECUTE.

****Compute unstandardized factor scores for p from bifactor model

COMPUTE p_fs_bi3=SUM(comp1_bi_p_weight3,q16_bi_p_weight3,comp2_bi_p_weight3,q37_bi_p_weight3,
    q95_bi_p_weight3,q3_bi_p_weight3,comp3_bi_p_weight3,q68_bi_p_weight3,q26_bi_p_weight3,q90_bi_p_weight3,
    q94_bi_p_weight3,comp4_bi_p_weight3,q86_bi_p_weight3,q43_bi_p_weight3,q67_bi_p_weight3,q87_bi_p_weight3,
    q27_bi_p_weight3,comp5_bi_p_weight3,q89_bi_p_weight3,q19_bi_p_weight3,q39_bi_p_weight3,q34_bi_p_weight3,
    q88_bi_p_weight3,q7_bi_p_weight3,q109_bi_p_weight3,q50_bi_p_weight3,q112_bi_p_weight3,q32_bi_p_weight3,q52_bi_p_weight3,
    q45_bi_p_weight3,q31_bi_p_weight3,q35_bi_p_weight3,q30_bi_p_weight3,q29_bi_p_weight3,q12_bi_p_weight3,
    comp6_bi_p_weight3,q17_bi_p_weight3,q13_bi_p_weight3,q62_bi_p_weight3,
    q46_bi_p_weight3,q4_bi_p_weight3,q9_bi_p_weight3,q61_bi_p_weight3,q66_bi_p_weight3,q1_bi_p_weight3,
    q36_bi_p_weight3,q64_bi_p_weight3, q56c_bi_p_weight3,q56f_bi_p_weight3,q56g_bi_p_weight3,q56b_bi_p_weight3,
    q56a_bi_p_weight3,q51_bi_p_weight3,q56d_bi_p_weight3,q56e_bi_p_weight3,
    q111_bi_p_weight3,q42_bi_p_weight3,q75_bi_p_weight3,q65_bi_p_weight3,
    q102_bi_p_weight3).
EXECUTE.
