%% aerodynamics look-up
function [ReList,wingPolar,wingPolarp20,wingPolarm20,tailPolarp30,A0,A0_tail,s0, s0_tail,s0_LE,s0_TE,Iy_LE,Ix_LE,Iy_TE,Ix_TE]=init
ReList=[30000,60000,100000,150000,200000,250000,300000,350000,400000,500000,600000];
% here's the polars
wingPolarp20=[0.000	1.073	0.072	-0.191];
wingPolarm20=[0.000	-0.794	0.045	0.064];
tailPolarp30=[0.000	1.163	0.058	-0.147];
wingPolar{1}=[...
-8	-0.2238	0.10231	-0.006
-7.5	-0.1704	0.09188	-0.0082
-7	-0.0737	0.07742	-0.0214
-6.5	0.0126	0.06464	-0.039
-6	0.0269	0.06029	-0.0399
-5.5	0.0389	0.0563	-0.0403
-4.5	-0.1318	0.06385	0.0004
-3	-0.3312	0.05348	-0.0241
-2.5	-0.2253	0.04296	-0.0458
-2	-0.1565	0.03827	-0.0502
-1.5	-0.0995	0.03495	-0.0482
-0.5	0.0025	0.03101	-0.0481
0	0.0468	0.0325	-0.0494
0.5	0.0886	0.03423	-0.0505
1	0.1285	0.03618	-0.0515
1.5	0.1663	0.03839	-0.0525
2	0.2023	0.04086	-0.0535
2.5	0.2363	0.04361	-0.0546
3	0.2684	0.04666	-0.0558
3.5	0.2986	0.05004	-0.0571
4	0.3269	0.05375	-0.0585
4.5	0.3536	0.05779	-0.0601
5	0.3789	0.06215	-0.0618
5.5	0.4029	0.06683	-0.0636
6	0.426	0.0718	-0.0655
6.5	0.4482	0.07705	-0.0674
7	0.4697	0.08256	-0.0695
7.5	0.4907	0.08832	-0.0717
8	0.511	0.09432	-0.0739
8.5	0.5561	0.1032	-0.0819
9	0.6064	0.11266	-0.0903
9.5	0.639	0.12029	-0.0944
10	0.671	0.12804	-0.098
10.5	0.7028	0.13626	-0.1013
11	0.7289	0.1439	-0.1036
11.5	0.7484	0.1506	-0.1049
12	0.7678	0.15752	-0.1064
12.5	0.7928	0.16566	-0.1086
13	0.8189	0.17396	-0.1108
13.5	0.8306	0.17924	-0.1118
14	0.8575	0.18798	-0.1141
14.5	0.8748	0.19473	-0.1158
15	0.8899	0.20056	-0.1176
15.5	0.9081	0.20723	-0.1197
16	0.9275	0.2144	-0.122
];
wingPolar{2}=[...
-8	-0.3169	0.1018	-0.0273
-7.5	-0.3475	0.09967	-0.0243
-7	-0.3423	0.0944	-0.019
-6.5	-0.3863	0.09335	-0.0103
-6	-0.396	0.0895	-0.003
-5.5	-0.4168	0.0864	0.0037
-5	-0.4356	0.08291	0.0095
-4.5	-0.4408	0.07917	0.0167
-3.5	-0.4399	0.07101	0.0352
-2.5	-0.3471	0.05135	-0.0094
-2	-0.1176	0.03134	-0.0607
-1.5	-0.0607	0.02936	-0.0615
-0.5	0.0183	0.02685	-0.054
0	0.0626	0.02858	-0.0556
0.5	0.1041	0.03054	-0.0568
1	0.1431	0.03273	-0.058
1.5	0.18	0.03518	-0.0592
2	0.2148	0.03792	-0.0603
2.5	0.2476	0.04096	-0.0616
3	0.2784	0.04432	-0.0629
3.5	0.3075	0.04803	-0.0644
4	0.3349	0.05207	-0.066
4.5	0.4132	0.05884	-0.0779
5	0.4812	0.06506	-0.0865
5.5	0.5281	0.0703	-0.0903
6	0.5809	0.07609	-0.0946
6.5	0.6192	0.08157	-0.0965
7	0.6426	0.08681	-0.0965
7.5	0.6888	0.09368	-0.0995
8	0.6949	0.09867	-0.0977
8.5	0.7281	0.10518	-0.0989
9	0.7415	0.11115	-0.0983
9.5	0.7705	0.1178	-0.0989
10	0.7806	0.124	-0.0983
10.5	0.8109	0.13124	-0.0989
11	0.816	0.13751	-0.0984
11.5	0.8373	0.1445	-0.0986
12	0.8604	0.15251	-0.0991
12.5	0.8631	0.15852	-0.0996
13	0.8755	0.16537	-0.1005
13.5	0.8922	0.1727	-0.1016
14	0.919	0.1813	-0.1025
14.5	0.9203	0.18631	-0.1047
15	0.9332	0.19353	-0.1069
15.5	0.947	0.20002	-0.1092
16	0.9752	0.20983	-0.1104

];
wingPolar{3}=[...
-7	-0.4211	0.09185	-0.0238
-4.5	-0.4385	0.06447	-0.027
-2.5	-0.1759	0.03089	-0.0586
-2	-0.1214	0.02866	-0.0596
-1.5	0.0186	0.02626	-0.0749
-1	0.0854	0.02403	-0.0739
-0.5	0.1788	0.02535	-0.0829
0	0.257	0.02648	-0.0886
0.5	0.3329	0.0275	-0.0933
1	0.4084	0.02832	-0.0973
2	0.5611	0.029	-0.1035
2.5	0.6414	0.02863	-0.1061
3	0.7171	0.02791	-0.1071
4	0.8596	0.02585	-0.1063
4.5	0.9247	0.02485	-0.1049
5	0.9856	0.02395	-0.103
5.5	1.0484	0.02282	-0.1015
6	1.1064	0.0222	-0.0997
6.5	1.1555	0.02243	-0.0972
7	1.2072	0.02287	-0.0952
7.5	1.2517	0.02382	-0.0925
8	1.2959	0.0249	-0.0899
8.5	1.3357	0.02631	-0.0869
9	1.3719	0.02789	-0.0835
9.5	1.402	0.02987	-0.0792
10	1.4309	0.03224	-0.0751
10.5	1.4539	0.03509	-0.0703
11	1.4799	0.03829	-0.0662
11.5	1.5062	0.04216	-0.0624
12	1.5344	0.04686	-0.0593
12.5	1.5478	0.05209	-0.0546
13	1.5641	0.0582	-0.051
13.5	1.5198	0.0645	-0.0412
14	1.4677	0.07284	-0.0348
14.5	1.4081	0.08371	-0.0324
15	1.3501	0.09723	-0.0339

];
wingPolar{4}=[...
-7.5	-0.4235	0.09052	-0.0299
-7	-0.4442	0.08833	-0.0225
-6.5	-0.4655	0.08417	-0.0213
-6	-0.4752	0.07347	-0.0343
-5.5	-0.4146	0.06387	-0.0467
-4.5	-0.2894	0.04818	-0.067
-4	-0.2219	0.04289	-0.0734
-3	-0.0035	0.02402	-0.0926
-2.5	0.0765	0.02176	-0.0971
-1.5	0.2088	0.01798	-0.0988
-1	0.2867	0.0174	-0.1012
-0.5	0.3747	0.01709	-0.1076
0	0.4411	0.01709	-0.1095
1	0.5794	0.01634	-0.1122
1.5	0.6347	0.01638	-0.1112
2	0.694	0.01622	-0.1106
2.5	0.7525	0.01608	-0.1098
3	0.807	0.01615	-0.1085
3.5	0.8617	0.01619	-0.1072
4	0.9168	0.0161	-0.1059
4.5	0.9673	0.01598	-0.1039
5	1.0204	0.01575	-0.1023
5.5	1.071	0.0158	-0.1004
6	1.1204	0.0162	-0.0985
6.5	1.1651	0.01693	-0.096
7	1.2095	0.01783	-0.0937
7.5	1.2505	0.01888	-0.091
8	1.2904	0.02007	-0.0883
8.5	1.3276	0.02137	-0.0853
9	1.3607	0.02281	-0.0817
9.5	1.3884	0.0245	-0.0775
10	1.4058	0.02645	-0.0718
10.5	1.4179	0.02879	-0.0655
11	1.431	0.03149	-0.06
11.5	1.4502	0.03456	-0.0557
12	1.4787	0.03798	-0.0529
12.5	1.501	0.04158	-0.0494
13	1.5252	0.04567	-0.0465
13.5	1.5451	0.05047	-0.0435
14	1.5753	0.05622	-0.042
14.5	1.5408	0.0622	-0.0356
15	1.5044	0.07013	-0.0318
15.5	1.4587	0.08006	-0.0305
16	1.3981	0.09303	-0.0326

];
wingPolar{5}=[...
-7.5	-0.4514	0.08993	-0.029
-6	-0.3276	0.05501	-0.0748
-3.5	0.0501	0.01936	-0.1076
-3	0.1218	0.01714	-0.1106
-2.5	0.191	0.01396	-0.1141
-2	0.2499	0.01316	-0.1135
-1.5	0.3008	0.01292	-0.1107
-1	0.339	0.01256	-0.1044
-0.5	0.4315	0.01206	-0.1111
0	0.4868	0.01208	-0.1112
0.5	0.5431	0.01211	-0.1109
1	0.5999	0.01217	-0.1105
1.5	0.6548	0.01229	-0.1098
2	0.7091	0.01248	-0.1089
2.5	0.7655	0.01262	-0.1083
3	0.8188	0.0128	-0.1073
3.5	0.8717	0.01288	-0.1061
4	0.9242	0.01283	-0.1048
4.5	0.9752	0.01279	-0.1032
5	1.0254	0.01292	-0.1016
5.5	1.074	0.01331	-0.0998
6	1.1195	0.01394	-0.0976
6.5	1.1635	0.01474	-0.0953
7	1.2056	0.01565	-0.0928
7.5	1.2467	0.01672	-0.0903
8	1.2859	0.01784	-0.0877
8.5	1.3223	0.01903	-0.0846
9	1.3548	0.02033	-0.081
9.5	1.3809	0.02182	-0.0766
10	1.3965	0.0236	-0.0706
10.5	1.4075	0.02598	-0.0647
11	1.4182	0.02896	-0.0595
11.5	1.4314	0.03202	-0.0551
12	1.4471	0.03529	-0.0513
12.5	1.4646	0.03865	-0.048
13	1.4837	0.0422	-0.045
13.5	1.5046	0.04601	-0.0422
14	1.5284	0.05012	-0.0399
14.5	1.53	0.05483	-0.0368
15	1.5402	0.06053	-0.0345
15.5	1.5203	0.06756	-0.0322
16	1.5208	0.07337	-0.031

];

wingPolar{6}=[...
-7.5	-0.4075	0.08613	-0.0302
-4	0.0462	0.01688	-0.1146
-3.5	0.1096	0.01511	-0.1157
-3	0.1691	0.01358	-0.1161
-2.5	0.2225	0.01112	-0.1165
-2	0.2754	0.0107	-0.115
-1.5	0.3224	0.01052	-0.1121
-1	0.3647	0.01036	-0.1077
-0.5	0.4257	0.01001	-0.1067
0	0.4947	0.00997	-0.1097
0.5	0.5481	0.01014	-0.1093
1	0.6046	0.01028	-0.1091
1.5	0.6585	0.0105	-0.1085
2	0.7146	0.0107	-0.1081
2.5	0.7676	0.0109	-0.1072
3	0.8208	0.01103	-0.1063
3.5	0.8738	0.01104	-0.1053
4	0.9262	0.01104	-0.1042
4.5	0.9771	0.01117	-0.1028
5	1.0264	0.01148	-0.1012
5.5	1.0739	0.012	-0.0994
6	1.1192	0.01266	-0.0973
6.5	1.1626	0.01348	-0.0951
7	1.2049	0.01439	-0.0928
7.5	1.2461	0.0154	-0.0904
8	1.2858	0.01645	-0.0878
8.5	1.323	0.01759	-0.0849
9	1.3561	0.01886	-0.0815
9.5	1.3833	0.02021	-0.0772
10	1.4014	0.02198	-0.0717
10.5	1.4142	0.02438	-0.0662
11	1.4262	0.02715	-0.0612
11.5	1.4357	0.03045	-0.0567
12	1.4482	0.03373	-0.053
12.5	1.4602	0.03728	-0.0496
13	1.473	0.04099	-0.0465
13.5	1.4882	0.04474	-0.0436
14	1.501	0.04857	-0.0413
14.5	1.5099	0.05314	-0.0388
15	1.529	0.0574	-0.0364
15.5	1.5206	0.06335	-0.0348
16	1.5422	0.06815	-0.0328


];

wingPolar{7}=[...
-8	-0.3294	0.07504	-0.0565
-7	-0.3231	0.06258	-0.0718
-3.5	0.1191	0.0135	-0.1152
-3	0.1731	0.01212	-0.1148
-2.5	0.2245	0.00998	-0.1148
-2	0.2773	0.00959	-0.1137
-1.5	0.3282	0.00926	-0.1118
-1	0.3765	0.00917	-0.1091
-0.5	0.4148	0.00896	-0.1036
0	0.4943	0.00886	-0.1083
0.5	0.55	0.00901	-0.1083
1	0.6047	0.00921	-0.1081
1.5	0.6601	0.00942	-0.1077
2	0.7146	0.00963	-0.1072
2.5	0.7693	0.00977	-0.1067
3	0.8223	0.00982	-0.1059
3.5	0.8752	0.00985	-0.1049
4	0.9275	0.00996	-0.1039
4.5	0.9778	0.01019	-0.1026
5	1.027	0.01061	-0.1011
5.5	1.0742	0.01115	-0.0994
6	1.1195	0.01185	-0.0974
6.5	1.1634	0.01265	-0.0953
7	1.2062	0.01356	-0.0931
7.5	1.248	0.01453	-0.0909
8	1.2884	0.01554	-0.0885
8.5	1.3272	0.01656	-0.0859
9	1.3623	0.01772	-0.0828
9.5	1.391	0.01905	-0.0787
10	1.413	0.0207	-0.0738
10.5	1.4279	0.02299	-0.0684
11	1.4397	0.02572	-0.0634
11.5	1.4484	0.02896	-0.0587
12	1.4595	0.03227	-0.0549
12.5	1.4686	0.03597	-0.0515
13	1.4769	0.03994	-0.0484
13.5	1.4866	0.04386	-0.0459
14	1.4946	0.04816	-0.0435
14.5	1.5043	0.05244	-0.0409
15	1.5074	0.05735	-0.0395
16	1.5123	0.06812	-0.0367

];

wingPolar{8}=[...
-5	-0.0464	0.01896	-0.1151
-4.5	0.0089	0.01521	-0.1146
-4	0.0643	0.01372	-0.1142
-3.5	0.118	0.01219	-0.1137
-3	0.1728	0.01124	-0.1134
-2.5	0.2251	0.0093	-0.1135
-2	0.2785	0.00885	-0.1128
-1.5	0.3315	0.00857	-0.1117
-1	0.3805	0.00844	-0.1094
-0.5	0.4226	0.00829	-0.1051
0	0.4938	0.00812	-0.1073
0.5	0.5502	0.0083	-0.1076
1	0.6048	0.00848	-0.1074
1.5	0.6604	0.00868	-0.1072
2	0.715	0.00884	-0.1068
2.5	0.7695	0.00893	-0.1063
3	0.8228	0.00897	-0.1055
3.5	0.8758	0.00904	-0.1047
4	0.9277	0.00922	-0.1037
4.5	0.9783	0.00955	-0.1025
5	1.0275	0.00999	-0.1011
5.5	1.0749	0.01056	-0.0995
6	1.1207	0.01124	-0.0977
6.5	1.1649	0.01206	-0.0957
7	1.2088	0.01292	-0.0938
7.5	1.2514	0.01384	-0.0917
8	1.2923	0.01484	-0.0894
8.5	1.3326	0.0158	-0.087
9	1.3693	0.0169	-0.0842
9.5	1.4003	0.01816	-0.0804
10	1.4237	0.01979	-0.0756
10.5	1.44	0.02198	-0.0703
11	1.4532	0.02455	-0.0653
11.5	1.4602	0.0278	-0.0604
12	1.4717	0.03101	-0.0566
12.5	1.4783	0.03485	-0.053
13	1.4901	0.03845	-0.0505
13.5	1.4964	0.04272	-0.0479
14	1.4993	0.04738	-0.0452
14.5	1.505	0.05198	-0.0437
15	1.508	0.05693	-0.0415
15.5	1.5094	0.06239	-0.0408
16	1.5138	0.06744	-0.039

];
wingPolar{9}=[...
-7.5	-0.3438	0.07114	-0.0566
-7	-0.2891	0.04786	-0.0929
-4.5	0.0082	0.01402	-0.1136
-4	0.0625	0.01256	-0.113
-3.5	0.1171	0.01157	-0.1126
-3	0.1718	0.01063	-0.1123
-2.5	0.2249	0.00894	-0.1125
-2	0.2784	0.00835	-0.112
-1.5	0.3321	0.00813	-0.1112
-1	0.3827	0.00793	-0.1095
-0.5	0.4293	0.00782	-0.1065
0	0.4791	0.00767	-0.1039
0.5	0.5496	0.00777	-0.1071
1	0.6048	0.00794	-0.107
1.5	0.6602	0.0081	-0.1068
2	0.7154	0.00824	-0.1066
2.5	0.7695	0.00828	-0.106
3	0.8231	0.00834	-0.1054
3.5	0.8762	0.00846	-0.1046
4	0.9281	0.0087	-0.1037
4.5	0.9788	0.00906	-0.1026
5	1.0282	0.00953	-0.1013
5.5	1.0758	0.01012	-0.0997
6	1.1221	0.0108	-0.098
6.5	1.1671	0.0116	-0.0962
7	1.2116	0.01244	-0.0944
7.5	1.255	0.01333	-0.0925
8	1.297	0.01426	-0.0904
8.5	1.3382	0.01519	-0.0882
9	1.3759	0.01626	-0.0855
9.5	1.4089	0.0175	-0.0821
10	1.433	0.01908	-0.0773
10.5	1.4505	0.02118	-0.072
11	1.4651	0.0236	-0.067
11.5	1.478	0.02634	-0.0625
12	1.4885	0.0295	-0.0584
12.5	1.4991	0.03292	-0.0551
13	1.5047	0.03701	-0.052
13.5	1.5043	0.04182	-0.049
14	1.514	0.04592	-0.0474
14.5	1.512	0.05119	-0.0451
15	1.5155	0.05631	-0.0441
15.5	1.5137	0.0619	-0.0424
16	1.5127	0.0681	-0.0426

];
wingPolar{10}=[...
-7	-0.2743	0.04173	-0.0994
-5.5	-0.1029	0.0177	-0.1132
-5	-0.0473	0.01464	-0.1124
-4.5	0.0053	0.01252	-0.1117
-4	0.0593	0.0114	-0.1113
-3.5	0.1143	0.01063	-0.1109
-3	0.1692	0.00976	-0.1108
-2.5	0.2236	0.00856	-0.1109
-2	0.2772	0.0077	-0.1109
-1.5	0.3321	0.00753	-0.1105
-1	0.3856	0.00729	-0.1096
-0.5	0.4362	0.00722	-0.1078
0	0.4814	0.00709	-0.1046
0.5	0.548	0.007	-0.1063
1	0.6043	0.00714	-0.1065
1.5	0.6604	0.00726	-0.1065
2	0.7152	0.00731	-0.1062
2.5	0.77	0.00739	-0.1058
3	0.8239	0.0075	-0.1053
3.5	0.8771	0.0077	-0.1046
4	0.929	0.00801	-0.1038
4.5	0.9799	0.00842	-0.1028
5	1.0299	0.00889	-0.1017
5.5	1.0781	0.00949	-0.1003
6	1.1253	0.01018	-0.0988
6.5	1.1713	0.01095	-0.0972
7	1.2169	0.01175	-0.0957
7.5	1.2617	0.01257	-0.094
8	1.3059	0.01339	-0.0922
8.5	1.3475	0.01433	-0.0901
9	1.3877	0.0153	-0.0878
9.5	1.4225	0.01654	-0.0847
10	1.4471	0.01811	-0.0799
10.5	1.4676	0.01998	-0.0748
11	1.486	0.02205	-0.07
11.5	1.5005	0.02455	-0.0653
12	1.5142	0.02731	-0.0612
12.5	1.5279	0.03032	-0.0578
13	1.5347	0.03416	-0.0546
13.5	1.5334	0.03897	-0.0517
14	1.5449	0.04286	-0.0499
14.5	1.5352	0.04888	-0.0476
15	1.5413	0.05375	-0.0468
15.5	1.5339	0.06022	-0.0459
16	1.53	0.06674	-0.0458

];
wingPolar{11}=[...
-5.5	-0.1044	0.0147	-0.1117
-5	-0.0515	0.0127	-0.1109
-4.5	0.0027	0.01161	-0.1105
-4	0.0573	0.0107	-0.1102
-3	0.1675	0.00922	-0.1099
-2.5	0.2227	0.0083	-0.11
-2	0.2766	0.00733	-0.1102
-1.5	0.3318	0.00713	-0.11
-1	0.3867	0.00696	-0.1095
-0.5	0.4389	0.00677	-0.1083
0	0.4885	0.00668	-0.1063
0.5	0.5366	0.00653	-0.1037
1	0.604	0.00656	-0.1062
1.5	0.6602	0.00664	-0.1062
2	0.7156	0.00669	-0.1061
2.5	0.7705	0.00679	-0.1058
3	0.8247	0.00695	-0.1053
3.5	0.8779	0.00721	-0.1047
4	0.9303	0.00754	-0.104
4.5	0.9816	0.00796	-0.1032
5	1.0317	0.00847	-0.1021
5.5	1.0806	0.00906	-0.1009
6	1.1283	0.00975	-0.0996
6.5	1.1755	0.01047	-0.0982
7	1.2221	0.01122	-0.0968
7.5	1.2677	0.01201	-0.0953
8	1.3128	0.01279	-0.0936
8.5	1.3557	0.01367	-0.0917
9	1.3961	0.01467	-0.0894
9.5	1.4317	0.01589	-0.0865
10	1.4591	0.01736	-0.0821
10.5	1.481	0.01909	-0.0771
11	1.5021	0.02094	-0.0724
11.5	1.521	0.02303	-0.0679
12	1.538	0.02543	-0.0639
12.5	1.5512	0.02832	-0.0601
13	1.5572	0.03207	-0.0565
13.5	1.5713	0.03539	-0.0542
14	1.5726	0.04009	-0.0517
14.5	1.5759	0.04481	-0.0499
15	1.5693	0.0508	-0.0483
15.5	1.5694	0.05647	-0.0477
16	1.5498	0.06464	-0.0473

];

airfoil = [...
  0.99991827     0.00088256
  0.99239612     0.00372417
  0.98156535     0.00691616
  0.97059460     0.00959046
  0.95957466     0.01205572
  0.94853628     0.01443719
  0.93748854     0.01677492
  0.92643480     0.01908410
  0.91537940     0.02138534
  0.90432443     0.02368861
  0.89326857     0.02598764
  0.88220902     0.02826881
  0.87114419     0.03052428
  0.86007493     0.03275788
  0.84899982     0.03496231
  0.83791622     0.03712359
  0.82682490     0.03924495
  0.81572749     0.04133424
  0.80462328     0.04338702
  0.79351176     0.04539990
  0.78239401     0.04737813
  0.77127065     0.04932421
  0.76014125     0.05123683
  0.74900790     0.05312119
  0.73786937     0.05499303
  0.72674267     0.05687245
  0.72000000     0.05803274
  0.71558714     0.05879211
  0.70456078     0.06073576
  0.69020760     0.06329651
  0.66935928     0.06690547
  0.64355338     0.07108051
  0.61581231     0.07512173
  0.58740137     0.07880793
  0.55883764     0.08206024
  0.53006340     0.08488644
  0.50123643     0.08733758
  0.47245118     0.08934950
  0.44373367     0.09092234
  0.41502894     0.09203455
  0.38639368     0.09266886
  0.35782287     0.09282508
  0.32933849     0.09247330
  0.30100774     0.09159539
  0.27286158     0.09017246
  0.24498189     0.08816427
  0.21751153     0.08553008
  0.19058666     0.08221354
  0.16432978     0.07813746
  0.13888591     0.07322849
  0.11445833     0.06743957
  0.09138306     0.06078041
  0.07025512     0.05340018
  0.05195427     0.04570434
  0.03725333     0.03828816
  0.02615023     0.03156807
  0.01798908     0.02563645
  0.01200709     0.02038662
  0.00761477     0.01567061
  0.00441337     0.01137293
  0.00216352     0.00739939
  0.00072642     0.00369208
  0.00003409     0.00023073
  0.00008807    -0.00302417
  0.00115976    -0.00624116
  0.00350412    -0.00909020
  0.00694098    -0.01147141
  0.01137268    -0.01357069
  0.01700203    -0.01550975
  0.02430623    -0.01735812
  0.03413767    -0.01918157
  0.04783529    -0.02102630
  0.06673654    -0.02281758
  0.09059551    -0.02430222
  0.11758053    -0.02529532
  0.14610326    -0.02584887
  0.17532475    -0.02606096
  0.20488466    -0.02600007
  0.23461316    -0.02570678
  0.26440235    -0.02520738
  0.29421282    -0.02451659
  0.32400889    -0.02362421
  0.35372912    -0.02252543
  0.38349605    -0.02119682
  0.41349344    -0.01965946
  0.44369834    -0.01801663
  0.47405359    -0.01630271
  0.50438236    -0.01455974
  0.53464892    -0.01282857
  0.56472804    -0.01108889
  0.59461297    -0.00950939
  0.62429556    -0.00807239
  0.65295448    -0.00687224
  0.67823259    -0.00592599
  0.69714806    -0.00528184
  0.70967055    -0.00483484
  0.72078924    -0.00438006
  0.73211349    -0.00384362
  0.74337826    -0.00324842
  0.75465675    -0.00262843
  0.76593234    -0.00202853
  0.77721124    -0.00147271
  0.78849189    -0.00095981
  0.79977435    -0.00048663
  0.81105891    -0.00006692
  0.82234582     0.00028408
  0.83363455     0.00057068
  0.84492451     0.00080378
  0.85621549     0.00097996
  0.86750726     0.00109615
  0.87879940     0.00116768
  0.89009170     0.00120539
  0.90138406     0.00119857
  0.91267626     0.00113855
  0.92396805     0.00102504
  0.93525910     0.00085362
  0.94654906     0.00062090
  0.95783739     0.00031949
  0.96912357    -0.00005384
  0.98040833    -0.00046804
  0.99169424    -0.00087171
  0.99974397    -0.00111748    
  0.99991827     0.00088256
];

LE_start=52;
LE_end=74;
leadingEdge=[airfoil(LE_start:LE_end,:);airfoil(LE_start,:)];
leadingEdge(LE_end-LE_start+1,:)=airfoil(LE_end-1,:)...
    +(airfoil(LE_end,:)-airfoil(LE_end-1,:))...
    *(1-(airfoil(LE_end,1)-(airfoil(LE_start,1)))/(airfoil(LE_end,1)-airfoil(LE_end-1,1)));

TE_start=110;
TE_end=13;
trailingEdge =[airfoil(TE_start:length(airfoil),:);...
    airfoil(1:TE_end,:);airfoil(TE_start,:)];
trailingEdge(length(trailingEdge)-1,:)=airfoil(TE_end-1,:)...
    +(airfoil(TE_end,:)-airfoil(TE_end-1,:))...
    *(1-(airfoil(TE_end,1)-(airfoil(TE_start,1)))/(airfoil(TE_end,1)-airfoil(TE_end-1,1)));

% plot(airfoil(:,1),airfoil(:,2)), axis equal; hold on;
% plot(leadingEdge(:,1),leadingEdge(:,2),'r'), axis equal; hold on;
% plot(trailingEdge(:,1),trailingEdge(:,2),'r'), axis equal; hold on;

airfoil_tail = [...
1.000000    0.000000
     0.993212    0.000476
     0.980807    0.001348
     0.965749    0.002337
     0.948812    0.003479
     0.930913    0.004785
     0.912504    0.006222
     0.893820    0.007746
     0.875003    0.009326
     0.856178    0.010927
     0.837509    0.012496
     0.819120    0.014128
     0.800992    0.015826
     0.783625    0.017594
     0.766711    0.019735
     0.749257    0.022299
     0.731221    0.024933
     0.712835    0.027591
     0.694255    0.030335
     0.675672    0.033052
     0.657193    0.035681
     0.638771    0.038195
     0.620392    0.040628
     0.602311    0.042909
     0.584260    0.044957
     0.565823    0.046905
     0.547272    0.048861
     0.528965    0.050725
     0.510843    0.052400
     0.492706    0.053913
     0.474615    0.055284
     0.456618    0.056473
     0.438609    0.057476
     0.420554    0.058307
     0.402436    0.058980
     0.384358    0.059512
     0.366410    0.059855
     0.348413    0.059991
     0.330339    0.059971
     0.312327    0.059791
     0.294354    0.059416
     0.276311    0.058860
     0.258325    0.058164
     0.240500    0.057269
     0.222614    0.056161
     0.204830    0.054912
     0.187271    0.053446
     0.169754    0.051745
     0.152439    0.049871
     0.135466    0.047743
     0.118700    0.045376
     0.102621    0.042846
     0.087235    0.040012
     0.072890    0.037108
     0.060163    0.034152
     0.049137    0.031235
     0.039944    0.028518
     0.032385    0.025938
     0.026136    0.023572
     0.020972    0.021409
     0.016688    0.019303
     0.013074    0.017239
     0.009974    0.015271
     0.007314    0.013377
     0.005083    0.011465
     0.003287    0.009468
     0.001918    0.007384
     0.000944    0.005252
     0.000328    0.003119
     0.000035    0.001021
     0.000035   -0.001021
     0.000328   -0.003118
     0.000944   -0.005252
     0.001918   -0.007384
     0.003287   -0.009468
     0.005083   -0.011465
     0.007314   -0.013377
     0.009974   -0.015271
     0.013074   -0.017238
     0.016688   -0.019303
     0.020973   -0.021408
     0.026136   -0.023572
     0.032386   -0.025938
     0.039945   -0.028518
     0.049138   -0.031234
     0.060164   -0.034151
     0.072892   -0.037108
     0.087237   -0.040012
     0.102623   -0.042845
     0.118702   -0.045376
     0.135468   -0.047742
     0.152441   -0.049870
     0.169756   -0.051744
     0.187274   -0.053445
     0.204832   -0.054912
     0.222617   -0.056161
     0.240502   -0.057269
     0.258328   -0.058164
     0.276313   -0.058860
     0.294356   -0.059416
     0.312330   -0.059791
     0.330341   -0.059971
     0.348415   -0.059991
     0.366412   -0.059855
     0.384360   -0.059512
     0.402438   -0.058980
     0.420555   -0.058307
     0.438611   -0.057476
     0.456619   -0.056473
     0.474616   -0.055284
     0.492707   -0.053913
     0.510845   -0.052401
     0.528966   -0.050725
     0.547273   -0.048861
     0.565824   -0.046905
     0.584261   -0.044957
     0.602312   -0.042909
     0.620393   -0.040628
     0.638772   -0.038194
     0.657195   -0.035681
     0.675674   -0.033051
     0.694256   -0.030335
     0.712837   -0.027590
     0.731224   -0.024933
     0.749260   -0.022298
     0.766713   -0.019735
     0.783626   -0.017594
     0.800992   -0.015825
     0.819121   -0.014128
     0.837510   -0.012495
     0.856179   -0.010927
     0.875004   -0.009326
     0.893821   -0.007745
     0.912505   -0.006222
     0.930914   -0.004785
     0.948813   -0.003479
     0.965749   -0.002337
     0.980807   -0.001348
     0.993212   -0.000477
     1.000000    0.000000
];

%% get some geometrical properties
% airfoil area
A0=0;
% airfoil circumference
s0=0;
% centroid
c0=[0,0];
for i=1:size(airfoil,1)-1
    A0=A0+(airfoil(i,1)*airfoil(i+1,2)-airfoil(i+1,1)*airfoil(i,2));
    d=sqrt((airfoil(i+1,:)-airfoil(i,:))*(airfoil(i+1,:)-airfoil(i,:))');
    s0=s0+d;
    c0=c0+d*(airfoil(i+1,:)+airfoil(i,:))/2;
end

A0_tail=0;
s0_tail=0;
for i=1:size(airfoil_tail,1)-1
    A0_tail=A0_tail+(airfoil_tail(i,1)*airfoil_tail(i+1,2)-airfoil_tail(i+1,1)*airfoil_tail(i,2));
    d=sqrt((airfoil_tail(i+1,:)-airfoil_tail(i,:))*(airfoil_tail(i+1,:)-airfoil_tail(i,:))');
    s0_tail=s0_tail+d;
end

% % second order moment devided by thickness
% Ix=0; Iy=0;
% for i=1:size(airfoil,1)-1
%     d=sqrt((airfoil(i+1,:)-airfoil(i,:))*(airfoil(i+1,:)-airfoil(i,:))');
%     P=(airfoil(i+1,:)+airfoil(i,:))/2;
%     rd=(P-c0);
%     Ix=Ix+rd(1)^2*d;
%     Iy=Iy+rd(2)^2*d;
% end

% calculate the front section
c0_LE=[0,0];
A0_LE=0;
s0_LE=0;
for i=1:size(leadingEdge,1)-1
    A0_LE=A0_LE+(leadingEdge(i,1)*leadingEdge(i+1,2)-leadingEdge(i+1,1)*leadingEdge(i,2));
    d=sqrt((leadingEdge(i+1,:)-leadingEdge(i,:))*(leadingEdge(i+1,:)-leadingEdge(i,:))');
    s0_LE=s0_LE+d;
    c0_LE=c0_LE+d*(leadingEdge(i+1,:)+leadingEdge(i,:))/2;
end
Ix_LE=0; %area second moment of inertia divided by shell thickness
Iy_LE=0; %area second moment of inertia divided by shell thickness
for i=1:size(leadingEdge,1)-1
    d=sqrt((leadingEdge(i+1,:)-leadingEdge(i,:))*(leadingEdge(i+1,:)-leadingEdge(i,:))');
    x=c0_LE(1)-leadingEdge(i,1);
    y=c0_LE(2)-leadingEdge(i,2);
    Iy_LE=Iy_LE+x^2*d;
    Ix_LE=Ix_LE+y^2*d;
end

% calculate the trailing section
c0_TE=[0,0];
A0_TE=0;
s0_TE=0;
for i=1:size(trailingEdge,1)-1
    A0_TE=A0_TE+(trailingEdge(i,1)*trailingEdge(i+1,2)-trailingEdge(i+1,1)*trailingEdge(i,2));
    d=sqrt((trailingEdge(i+1,:)-trailingEdge(i,:))*(trailingEdge(i+1,:)-trailingEdge(i,:))');
    s0_TE=s0_TE+d;
    c0_TE=c0_TE+d*(trailingEdge(i+1,:)+trailingEdge(i,:))/2;
end
Ix_TE=0; %area second moment of inertia divided by shell thickness
Iy_TE=0; %area second moment of inertia divided by shell thickness
for i=1:size(trailingEdge,1)-1
    d=sqrt((trailingEdge(i+1,:)-trailingEdge(i,:))*(trailingEdge(i+1,:)-trailingEdge(i,:))');
    x=c0_TE(1)-trailingEdge(i,1);
    y=c0_TE(2)-trailingEdge(i,2);
    Iy_TE=Iy_TE+x^2*d;
    Ix_TE=Ix_TE+y^2*d;
end

%% calculation of deflected control surface airfoils
% r_ctr = 0.3;
% a_main = 20/180*pi;
% a_tail = 30/180*pi;
% ctr=[1-r_ctr; 0]; %center of rotation
% airfoil_p=[];
% for i=1:size(airfoil,1)-1
%     if (airfoil(i,1)>=1-r_ctr)
%         p_new=[cos(a_main) sin(a_main); -sin(a_main) cos(a_main)]...
%             *(airfoil(i,:)'-ctr)+ctr;
%         if (p_new(1)>=1-r_ctr)
%             airfoil_p=[airfoil_p; p_new'];
%         end
%     else
%         airfoil_p=[airfoil_p; airfoil(i,:)];
%     end
% end
% airfoil_m=[];
% for i=1:size(airfoil,1)-1
%     if (airfoil(i,1)>=1-r_ctr)
%         p_new=[cos(-a_main) sin(-a_main); -sin(-a_main) cos(-a_main)]...
%             *(airfoil(i,:)'-ctr)+ctr;
%         if (p_new(1)>=1-r_ctr)
%             airfoil_m=[airfoil_m; p_new'];
%         end
%     else
%         airfoil_m=[airfoil_m; airfoil(i,:)];
%     end
% end
% airfoil_tail_p=[];
% for i=1:size(airfoil_tail,1)
%     if (airfoil_tail(i,1)>=1-r_ctr)
%         p_new=[cos(a_tail) sin(a_tail); -sin(a_tail) cos(a_tail)]...
%             *(airfoil_tail(i,:)'-ctr)+ctr;
%         if (p_new(1)>=1-r_ctr)
%             airfoil_tail_p=[airfoil_tail_p; p_new'];
%         end
%     else
%         airfoil_tail_p=[airfoil_tail_p; airfoil_tail(i,:)];
%     end
% end
% airfoil_tail_m=[];
% for i=1:size(airfoil_tail,1)
%     if (airfoil_tail(i,1)>=1-r_ctr)
%         p_new=[cos(-a_tail) sin(-a_tail); -sin(-a_tail) cos(-a_tail)]...
%             *(airfoil_tail(i,:)'-ctr)+ctr;
%         if (p_new(1)>=1-r_ctr)
%             airfoil_tail_m=[airfoil_tail_m; p_new'];
%         end
%     else
%         airfoil_tail_m=[airfoil_tail_m; airfoil_tail(i,:)];
%     end
% end