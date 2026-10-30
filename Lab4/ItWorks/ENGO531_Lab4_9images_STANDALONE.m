function ENGO531_Lab4_9images_STANDALONE()
% ONE standalone 9-image course bundle adjustment. No previous scripts,
% MAT result files, or .tie file required at runtime.
% Starting EOP/XYZ seeds below are frozen numerical values recovered from
% YOUR previously converged 1-6 and 7-9 solutions (not generic estimates).
% REQUIRED: calibration(4).PHO and full_matlab_bundle.p on the MATLAB path.
% The output is one joint adjustment, not two separately reported networks.
clc;
numberOfBundleDistances = 2; % 2=previously converged; 1=assignment-style test
assert(ismember(numberOfBundleDistances,[1,2]));
IOP=[2834 2156.7 3995.9];
[f,d]=uigetfile({'*.PHO;*.pho','Image observations (*.PHO)'}, ...
    'Select your updated calibration(4).PHO');
if isequal(f,0),return;end
fid=fopen(fullfile(d,f),'r');assert(fid>=0,'Cannot open PHO file');
c=textscan(fid,'%s%f%f%f');fclose(fid);
rawName=string(c{1});rawIm=c{2};rawXY=[c{3} c{4}];
folder=uigetdir(d,'Select output folder for the 9-image adjustment');
if isequal(folder,0),return;end
prefix=fullfile(folder,'ENGO531_JOINT_1to9');
% Embedded INITIAL coordinates from original all-nine .tie file:
Tall=[...
1001.0000000000 -65.0000000000 2704.0000000000 0.0000000000;
1002.0000000000 -65.0000000000 563.7165640000 -1102.3223370000;
1003.0000000000 -65.0000000000 700.9311250000 -1863.0683980000;
1004.0000000000 -65.0000000000 1803.3251600000 0.8478290000;
1005.0000000000 -65.0000000000 900.4215290000 0.8564620000;
1006.0000000000 -65.0000000000 0.0000000000 -0.0000000000;
1007.0000000000 -65.0000000000 -0.0000000000 -2213.0000000000;
1008.0000000000 -65.0000000000 897.9939010000 -2213.1344150000;
1009.0000000000 -65.0000000000 1801.3161330000 -2212.7508760000;
1010.0000000000 -65.0000000000 2704.0000000000 -2213.0000000000;
1011.0000000000 -65.0000000000 1386.4350380000 -1101.2461870000;
1012.0000000000 -65.0000000000 1352.0000000000 -2213.0000000000;
1013.0000000000 -0.0000000000 0.0000000000 -2213.0000000000;
1014.0000000000 708.7520640000 0.0000000000 2.6461000000;
1015.0000000000 745.5766560000 0.0000000000 7.8144560000;
1016.0000000000 1726.5523260000 0.0000000000 4.9383480000;
1017.0000000000 2478.0000000000 0.0000000000 0.0000000000;
1018.0000000000 2478.0000000000 0.0000000000 -2213.0000000000;
1019.0000000000 -0.0000000000 0.0000000000 0.0000000000;
1020.0000000000 765.9514240000 0.0000000000 -1114.7284510000;
1021.0000000000 765.5716440000 0.0000000000 -1176.6161440000;
1022.0000000000 757.7855040000 0.0000000000 -2233.7547800000;
1023.0000000000 2543.0000000000 2704.0000000000 0.0000000000;
1024.0000000000 2543.0000000000 621.9407580000 -1108.9803550000;
1025.0000000000 2543.0000000000 1998.4026530000 -1861.3554380000;
1026.0000000000 2543.0000000000 2511.0246980000 -1971.2490940000;
1027.0000000000 2543.0000000000 192.8660860000 -1863.8521480000;
1028.0000000000 2543.0000000000 704.9958180000 -1971.5686410000;
1029.0000000000 2543.0000000000 1802.7892430000 1.1315190000;
1030.0000000000 2543.0000000000 899.0647450000 1.6299320000;
1031.0000000000 2543.0000000000 0.0000000000 0.0000000000;
1032.0000000000 2543.0000000000 0.0000000000 -2213.0000000000;
1033.0000000000 2543.0000000000 2704.0000000000 -2213.0000000000;
1034.0000000000 2543.0000000000 1440.1308660000 -1114.5403100000;
1035.0000000000 2543.0000000000 900.4394820000 -2211.6277320000;
1036.0000000000 2543.0000000000 1802.1904730000 -2211.1706990000
];
% Embedded independently converged six-image reference geometry:
A.tieTable=[...
1001.0000000000 -65.0000000000 563.7165640961 -1102.3223372732;
1002.0000000000 -65.0000000000 700.9311252388 -1863.0683975687;
1003.0000000000 -65.0000000000 1803.3251595795 0.8478292145;
1004.0000000000 -65.0000000000 0.0000000000 -0.0000000000;
1005.0000000000 -65.0000000000 -0.0000000000 -2213.0000000000;
1006.0000000000 -65.0000000000 1801.3161329631 -2212.7508764662;
1007.0000000000 -65.0000000000 1386.4350376629 -1101.2461867423;
1008.0000000000 -0.0000000000 0.0000000000 -2213.0000000000;
1009.0000000000 745.5766561931 0.0000000000 7.8144559607;
1010.0000000000 1726.5523259212 0.0000000000 4.9383475265;
1011.0000000000 2478.0000000000 0.0000000000 0.0000000000;
1012.0000000000 2478.0000000000 0.0000000000 -2213.0000000000;
1013.0000000000 -0.0000000000 0.0000000000 0.0000000000;
1014.0000000000 765.9514240420 0.0000000000 -1114.7284512827;
1015.0000000000 765.5716438649 0.0000000000 -1176.6161443693
];
A.tie.XYZ=[...
188.6699672003 803.7899057131 -765.3399272993;
265.2167182017 965.5163462451 -1458.2863451954;
-157.8327045928 1987.9184019077 267.3988572013;
128.2937964219 207.3975816327 239.6316851993;
424.6025143043 282.4211242027 -1792.8147273727;
139.8938961681 2063.9871698207 -1769.2215276057;
61.0255157549 1613.4914382516 -754.1015421206;
485.8932722248 310.1748433180 -1799.8465113764;
922.6373216420 326.5790035324 379.0290846815;
1875.9193816511 468.6554375146 520.5368837136;
2609.4018400157 591.5070347553 630.8160489197;
2917.1670747564 675.1519950188 -1429.6515184905;
184.7733036964 232.9752934847 265.9661464014;
1089.5244644646 380.4899839839 -658.1207033405;
1097.7119889312 382.7530417994 -716.4774424293
];
A.photoTable=[...
1001.0000000000 1.0000000000 4347.4700000000 2355.5400000000;
1002.0000000000 1.0000000000 4037.5500000000 3472.4500000000;
1003.0000000000 1.0000000000 2295.1200000000 384.4200000000;
1004.0000000000 1.0000000000 5399.7300000000 646.4800000000;
1005.0000000000 1.0000000000 4995.2300000000 3894.6700000000;
1006.0000000000 1.0000000000 2234.6300000000 4065.5000000000;
1007.0000000000 1.0000000000 3003.5900000000 2345.4400000000;
1001.0000000000 2.0000000000 3334.4300000000 2169.4000000000;
1002.0000000000 2.0000000000 3096.5200000000 3222.3800000000;
1003.0000000000 2.0000000000 1676.6900000000 537.5800000000;
1004.0000000000 2.0000000000 4506.9900000000 221.2900000000;
1005.0000000000 2.0000000000 4210.0600000000 4002.2800000000;
1006.0000000000 2.0000000000 1776.0300000000 3339.3000000000;
1007.0000000000 2.0000000000 2210.4500000000 2066.5700000000;
1008.0000000000 3.0000000000 2713.5600000000 4019.3100000000;
1009.0000000000 3.0000000000 3702.0800000000 484.9900000000;
1010.0000000000 3.0000000000 4452.4100000000 716.2700000000;
1011.0000000000 3.0000000000 4831.0800000000 838.0000000000;
1012.0000000000 3.0000000000 4618.4400000000 3023.0900000000;
1001.0000000000 3.0000000000 2104.4700000000 2141.4800000000;
1002.0000000000 3.0000000000 1999.5500000000 3125.4900000000;
1003.0000000000 3.0000000000 1266.0300000000 702.2700000000;
1004.0000000000 3.0000000000 2694.9600000000 219.3200000000;
1005.0000000000 3.0000000000 2648.2200000000 4051.9900000000;
1006.0000000000 3.0000000000 1367.9100000000 3094.8700000000;
1007.0000000000 3.0000000000 1533.4900000000 2006.5400000000;
1013.0000000000 4.0000000000 1142.1900000000 135.0700000000;
1008.0000000000 4.0000000000 1315.8500000000 3967.1800000000;
1014.0000000000 4.0000000000 2552.2500000000 2145.6800000000;
1015.0000000000 4.0000000000 2550.3200000000 2237.6600000000;
1009.0000000000 4.0000000000 2563.4900000000 351.3700000000;
1010.0000000000 4.0000000000 3862.0400000000 565.2800000000;
1011.0000000000 4.0000000000 4593.3000000000 686.0500000000;
1012.0000000000 4.0000000000 4335.7900000000 3194.5100000000;
1004.0000000000 4.0000000000 1007.9900000000 140.1600000000;
1005.0000000000 4.0000000000 1202.0600000000 3989.8200000000;
1013.0000000000 5.0000000000 428.1800000000 223.2300000000;
1008.0000000000 5.0000000000 666.6200000000 3908.1800000000;
1014.0000000000 5.0000000000 1937.2500000000 2245.0500000000;
1015.0000000000 5.0000000000 1937.4700000000 2345.4800000000;
1009.0000000000 5.0000000000 1880.5700000000 267.7300000000;
1010.0000000000 5.0000000000 3729.9900000000 348.9400000000;
1011.0000000000 5.0000000000 5099.8400000000 414.1200000000;
1012.0000000000 5.0000000000 4710.0800000000 3926.0100000000;
1013.0000000000 6.0000000000 1047.5000000000 555.5000000000;
1008.0000000000 6.0000000000 1142.5000000000 2891.5000000000;
1014.0000000000 6.0000000000 1573.5000000000 1900.5000000000;
1015.0000000000 6.0000000000 1574.5000000000 1977.5000000000;
1009.0000000000 6.0000000000 1521.5000000000 453.5000000000;
1010.0000000000 6.0000000000 2463.5000000000 273.5000000000;
1011.0000000000 6.0000000000 3645.4900000000 59.5100000000;
1012.0000000000 6.0000000000 3343.4900000000 4127.5100000000
];
A.eop.XYZc=[...
-2213.0000000000 1711.0000000000 -645.0000000000;
-1946.0988801089 -310.7705081319 -699.1441797459;
-840.5120298146 -1546.6945227281 -592.8403741331;
384.1585094593 -2013.2488216521 -389.2205098396;
1776.8586324388 -1794.7665873275 -203.7450838703;
3945.7569099522 -928.0178974851 170.0571195704
];
A.eop.wpk=[...
90.0000000000 -90.0000000000 -0.0000000000;
86.6851278550 -55.1520708066 1.0201691076;
86.5336861904 -33.9926559040 4.9358729362;
83.9271379113 -20.3406680329 6.6644535095;
81.9568436324 5.0196108052 10.9117710803;
67.1920047196 44.5899571589 26.4372625937
];
A.IOP=IOP;
% Embedded independently converged right-wall reference geometry:
B.tieTable=[...
1001.0000000000 2543.0000000000 2701.0000000000 0.0000000000;
1002.0000000000 2543.0000000000 621.2507346784 -1108.9803545457;
1003.0000000000 2543.0000000000 1996.1854904774 -1861.3554384774;
1004.0000000000 2543.0000000000 2508.2387978024 -1971.2490937156;
1005.0000000000 2543.0000000000 192.6521077706 -1863.8521480601;
1006.0000000000 2543.0000000000 704.2136479125 -1971.5686413516;
1007.0000000000 2543.0000000000 1800.7891072895 1.1315191189;
1008.0000000000 2543.0000000000 898.0672618237 1.6299316814;
1009.0000000000 2543.0000000000 0.0000000000 0.0000000000;
1010.0000000000 2543.0000000000 0.0000000000 -2213.0000000000;
1011.0000000000 2543.0000000000 899.4404739296 -2211.6277320178;
1012.0000000000 2543.0000000000 1800.1910013358 -2211.1706993310;
1013.0000000000 2543.0000000000 2701.0000000000 -2213.0000000000;
1014.0000000000 2543.0000000000 1438.5330872893 -1114.5403097312
];
B.tie.XYZ=[...
2556.3504929625 2755.6981415244 618.0610160218;
2213.4706029534 773.1869498584 -529.7027369544;
2202.8936553120 2195.0903278080 -1123.8276760551;
2229.8656894742 2712.3882308681 -1190.0136964241;
2059.8101064605 400.3456389943 -1247.3981078678;
2083.5565633505 918.0232318297 -1309.8385524514;
2480.5832483910 1863.5682241641 560.8562001603;
2404.6285223595 966.0367148680 502.3457538938;
2330.6429214858 69.9233395376 441.7918888583;
1987.5145551967 233.0791267083 -1580.4136906953;
2057.1953548865 1129.3283987746 -1519.1049080107;
2129.8525564749 2027.4941126408 -1458.5369008602;
2207.2151904037 2922.5176420515 -1398.4739471941;
2267.5455266306 1586.5739623413 -481.1560961682
];
B.eop.XYZc=[...
4368.0000000000 -97.0000000000 -495.0000000000;
4541.1445500746 1138.4644274252 -436.5391514716;
4584.7550064554 2451.9505181398 -347.8278555472
];
B.eop.wpk=[...
90.0000000000 60.0000000000 0.0000000000;
95.5930025119 80.0571682200 -1.4479545015;
98.2575292226 106.3122191534 1.1257470441
];
B.IOP=IOP;
B.phoTable=[...
1001.0000000000 7.0000000000 4920.5000000000 669.5000000000;
1002.0000000000 7.0000000000 2273.5000000000 2217.5000000000;
1003.0000000000 7.0000000000 4026.5000000000 2990.5000000000;
1004.0000000000 7.0000000000 4507.5000000000 3011.5100000000;
1005.0000000000 7.0000000000 1548.5000000000 3494.5100000000;
1006.0000000000 7.0000000000 2411.5000000000 3463.5000000000;
1007.0000000000 7.0000000000 3985.9800000000 541.5000000000;
1008.0000000000 7.0000000000 2724.9800000000 370.5000000000;
1009.0000000000 7.0000000000 943.5000000000 130.9700000000;
1010.0000000000 7.0000000000 1211.5000000000 4104.9800000000;
1011.0000000000 7.0000000000 2691.4900000000 3721.5000000000;
1012.0000000000 7.0000000000 3793.5000000000 3439.9800000000;
1013.0000000000 7.0000000000 4647.5000000000 3225.9700000000;
1014.0000000000 7.0000000000 3445.5000000000 2138.0000000000;
1001.0000000000 8.0000000000 5180.5000000000 514.5100000000;
1002.0000000000 8.0000000000 1455.5000000000 2294.0000000000;
1003.0000000000 8.0000000000 3782.5000000000 3405.5100000000;
1004.0000000000 8.0000000000 4559.5000000000 3536.5000000000;
1005.0000000000 8.0000000000 732.5000000000 3484.5000000000;
1006.0000000000 8.0000000000 1646.5000000000 3619.5000000000;
1007.0000000000 8.0000000000 3621.5000000000 440.5100000000;
1008.0000000000 8.0000000000 1931.9800000000 361.5100000000;
1009.0000000000 8.0000000000 110.5100000000 280.5000000000;
1010.0000000000 8.0000000000 423.5000000000 4009.5000000000;
1011.0000000000 8.0000000000 1986.5000000000 3950.9600000000;
1012.0000000000 8.0000000000 3447.0000000000 3897.5000000000;
1013.0000000000 8.0000000000 4812.5000000000 3850.5100000000;
1014.0000000000 8.0000000000 2912.5100000000 2305.4800000000;
1001.0000000000 9.0000000000 5045.0000000000 200.9700000000;
1002.0000000000 9.0000000000 1452.5000000000 2040.0000000000;
1003.0000000000 9.0000000000 3341.5000000000 3387.5000000000;
1004.0000000000 9.0000000000 4228.5000000000 3750.5000000000;
1005.0000000000 9.0000000000 1021.5000000000 2907.5000000000;
1006.0000000000 9.0000000000 1578.5000000000 3153.5000000000;
1007.0000000000 9.0000000000 3162.4900000000 331.5000000000;
1008.0000000000 9.0000000000 1765.5000000000 429.5000000000;
1009.0000000000 9.0000000000 690.5000000000 504.5000000000;
1010.0000000000 9.0000000000 854.5500000000 3259.5700000000;
1011.0000000000 9.0000000000 1811.9900000000 3500.9800000000;
1012.0000000000 9.0000000000 3012.5000000000 3808.5000000000;
1013.0000000000 9.0000000000 4560.5000000000 4209.5000000000;
1014.0000000000 9.0000000000 2526.5000000000 2162.5000000000
];
assert(size(Tall,2)==4 && size(A.tie.XYZ,1)==15 && size(B.tie.XYZ,1)==14);
% Recover identity from actual input 3D approximations, not numericIDs.
% These are shared across the original all-image tie file and each subset's
% saved tieTable; allow 10 mm for R4/R1 length 2701 vs 2704 discrepancy.
ixA=matchXYZ(A.tieTable(:,2:4),Tall(:,2:4),1.0,'images 1-6');
ixB=matchXYZ(B.tieTable(:,2:4),Tall(:,2:4),10,'images 7-9');
fprintf('Matched %d left/front points and %d right-wall points to 9-image initial geometry.\n',numel(ixA),numel(ixB));
assert(isempty(intersect(ixA,ixB)),'Unexpected common point match between distinct faces');
% For the original all-nine input .tie the IDs are assigned to the stable
% first-appearance order of point names in the raw PHO, after removing
% targets with fewer than 2 images. Reconstruct and verify by the X/Y/Z shape.
rawOK=ismember(rawIm,1:9);rawName=rawName(rawOK);rawIm=rawIm(rawOK);rawXY=rawXY(rawOK,:);
[uniqueNames,~,group]=unique(rawName,'stable');count=accumarray(group,1);
retained=uniqueNames(count>=2);
assert(numel(retained)==size(Tall,1),'All-image .tie point count does not match raw PHO filtered to >=2 images.');
expectedID=(1001:1000+numel(retained))';
assert(isequal(Tall(:,1),expectedID),'All-image .tie target IDs are not the expected sequential IDs.');
[isFound,groupFiltered]=ismember(rawName,retained);
P=[expectedID(groupFiltered(isFound)),rawIm(isFound),rawXY(isFound,:)];
assert(all(ismember(1:9,P(:,2))),'A camera has no retained measurements.');

% Validated collinearity convention from geometry test:
% M = (Rx(omega)*Ry(phi)*Rz(kappa))'  (world to camera)
% x = xp - c*qx/qz; y = yp + c*qy/qz; camera y-axis direction -1.
IOP=[2834 2156.7 3995.9];
if isfield(A,'IOP'),IOP=double(A.IOP(:)');end
assert(norm(IOP-double(B.IOP(:)'))<0.1,'The two independent adjustments use different IOPs');
% Explicit validation of the six-image MAT, BEFORE joining the networks.
PA=A.photoTable;[ok,kA]=ismember(PA(:,1),A.tieTable(:,1));assert(all(ok));
err=[];
fprintf('\nREFERENCE REPROJECTION, converged 1-6 solution:\n');
for i=1:6
 mask=PA(:,2)==i;
 predicted=project(A.tie.XYZ(kA(mask),:),[A.eop.XYZc(i,:) A.eop.wpk(i,:)],IOP);
 e=hypot(predicted(:,1)-PA(mask,3),predicted(:,2)-PA(mask,4));
 fprintf(' Image %d: %.3f px radial RMS (%d points)\n',i,sqrt(mean(e.^2)),numel(e));
 if i>=2,err=[err;e];end %#ok<AGROW>
end
refRMS=sqrt(mean(err.^2));
fprintf(' Verified images 2-6 reference RMS = %.3f px\n',refRMS);
assert(isfinite(refRMS)&&refRMS<15,'Reference validation failed; do not run joint bundle.');

% Resection: solve camera 6 in the right-wall coordinate frame from R targets
% observed in original Image 6 (not present in the 7-9 MAT's photoTable).
[ok,loc]=ismember(Tall(ixB,1),Tall(:,1));assert(all(ok));
[has,j]=ismember(P(:,1),Tall(ixB,1));
sel=P(:,2)==6 & has;
obs6=P(sel,3:4);idxB=j(sel);
assert(numel(idxB)>=8,'Cannot find enough image-6 R-wall observations');
Xr=B.tie.XYZ(idxB,:);
seeds=[3486 -1652 -495 90 30 0;3942 -933 170 67 45 26; ...
       4000 -1300 -300 90 45 0;4200 -500 -500 90 55 0; ...
       B.eop.XYZc(1,:) B.eop.wpk(1,:)];
opts=optimset('Display','off','MaxIter',4000,'MaxFunEvals',25000,'TolX',1e-8,'TolFun',1e-8);
best=Inf;bestP=[];
for ii=1:size(seeds,1)
 init=[seeds(ii,1:3)/1000,seeds(ii,4:6)/30];
 [v,~,flag]=fminsearch(@(v)poseCost(v,Xr,obs6,IOP),init,opts);
 p=[v(1:3)*1000 v(4:6)*30];
 pr=project(Xr,p,IOP);
 er=sqrt(mean(sum((pr-obs6).^2,2)));
 fprintf('Image 6 resection seed %d: RMS %.3f px; exit flag %d\n',ii,er,flag);
 if isfinite(er) && er<best,best=er;bestP=p;end
end
fprintf('Best Image 6 resection RMS: %.3f px\n',best);
if isempty(bestP)||best>35
 fprintf(2,'STOP: Image 6 resection not accurate enough to register networks. No bundle run.\n');
 return;
end

% Frame transform: worldLeft = R*worldRight + translation.
ML=rotW2C(A.eop.wpk(6,:)); MR=rotW2C(bestP(4:6));
R=ML'*MR; translation=A.eop.XYZc(6,:)'-R*bestP(1:3)';
TallNew=Tall;
% For points missing from the saved 1-6 solution, intersect camera rays
% from the ORIGINAL image observations using the CONVERGED 1-6 EOPs.
% This is essential: a rigid transform of estimated face coordinates is
% not precise enough for the extra L/F points in Images 1-3.
extraLF=setdiff((1:size(Tall,1))',[ixA;ixB]);
assert(all(startsWith(retained(extraLF),"L")|startsWith(retained(extraLF),"F")), ...
    'An unexpected unmatched point is not an L/F target.');
fprintf('\nTRIANGULATING %d missing L/F points from calibrated camera rays:\n',numel(extraLF));
triResults=zeros(numel(extraLF),4);
for j=1:numel(extraLF)
    idx=extraLF(j); pid=Tall(idx,1);
    mask=P(:,1)==pid & P(:,2)<=6;
    assert(nnz(mask)>=2,'Target %s has fewer than two 1-6 observations',retained(idx));
    obs=P(mask,:);
    [X,linearRank]=triangulateWeightedLinear(obs,A.eop.XYZc,A.eop.wpk,IOP);
    assert(linearRank==3 && all(isfinite(X)), ...
        'Triangulation failed for target %s. Inspect stereo geometry.',retained(idx));
    % Refine 3D point alone, with all six-image EOPs FIXED.
    options=optimset('Display','off','MaxIter',1500,'MaxFunEvals',5000, ...
        'TolX',1e-9,'TolFun',1e-9);
    X=fminsearch(@(x) pointCost(x,obs,A.eop.XYZc,A.eop.wpk,IOP),X,options);
    errs=pointErrors(X,obs,A.eop.XYZc,A.eop.wpk,IOP);
    rms=sqrt(mean(errs.^2));
    fprintf('  %-5s images %s: radial RMS %8.3f px, XYZ [%9.2f %9.2f %9.2f] mm\n', ...
        char(retained(idx)),mat2str(obs(:,2)'),rms,X(1),X(2),X(3));
    TallNew(idx,2:4)=X;
    triResults(j,:)=[pid,nnz(mask),rms,max(errs)];
end
TallNew(ixA,2:4)=A.tie.XYZ;
TallNew(ixB,2:4)=(R*B.tie.XYZ'+translation)';
if any(triResults(:,3)>15)
    warning('Some new points have >15px triangulation RMS. Inspect corresponding image labels.');
end
writetable(array2table(triResults,'VariableNames', ...
  {'ID','NumberOfRays','RadialRMS_px','MaxError_px'}),[prefix '_triangulation.csv']);
% Original source EOPs for image numbering; replace with converged poses.
notes=[1 -2213 1711 -645 -90 0 90;2 -2279 561 -495 -75 0 90; ...
3 -1619 -1540 -495 -45 0 90;4 -440 -2027 -495 -25 0 90; ...
5 1239 -2545 -495 0 0 90;6 3486 -1652 -495 30 0 90; ...
7 4368 -97 -495 60 0 90;8 4901 411 -495 90 0 90;9 4455 2686 -495 120 0 90];
E=[notes(:,1),99*ones(9,1),notes(:,2:4),notes(:,7),notes(:,5),notes(:,6)];
E(1:6,3:5)=A.eop.XYZc; E(1:6,6:8)=A.eop.wpk;
for i=7:9
 k=i-6;
 E(i,3:5)=(R*B.eop.XYZc(k,:)'+translation)';
 Mnew=rotW2C(B.eop.wpk(k,:))*R';
 E(i,6:8)=rot2WPK(Mnew);
end
% Preflight per camera against ORIGINAL measurements (raw PHO not changed)
[ok,k]=ismember(P(:,1),TallNew(:,1));assert(all(ok));
summary=zeros(9,3);
for i=1:9
 mask=P(:,2)==i; pred=project(TallNew(k(mask),2:4),E(i,3:8),IOP);
 e=hypot(pred(:,1)-P(mask,3),pred(:,2)-P(mask,4));
 summary(i,:)=[i,nnz(mask),sqrt(mean(e.^2))];
end
fprintf('\nPRE-BUNDLE REPROJECTION RMS, each image:\n');
disp(array2table(summary,'VariableNames',{'Image','Count','RadialRMS_px'}));
writetable(array2table(summary,'VariableNames',{'Image','Count','RadialRMS_px'}),[prefix '_preflight.csv']);
if max(summary(:,3))>100 || any(~isfinite(summary(:,3)))
 fprintf(2,'STOP: initial geometry is not consistent enough. Joint bundle not run.\n');
 return;
end
% Create NEW course files, original image observations preserved.
scaleA=find(retained=="F1",1);scaleB=find(retained=="F6",1);
assert(~isempty(scaleA)&&~isempty(scaleB),'Scale endpoints F1/F6 absent');
writeRows([prefix '.pho'],P,'%d %d %.4f %.4f\n');
writeRows([prefix '.tie'],TallNew,'%d %.7f %.7f %.7f\n');
writeRows([prefix '.ext'],E,'%d %d %.7f %.7f %.7f %.8f %.8f %.8f\n');
writeRows([prefix '.int'],[99 IOP -1],'%d %.5f %.5f %.5f %d\n');
distanceRows=[expectedID(scaleA),expectedID(scaleB),2478];
if numberOfBundleDistances==2
    r4=find(retained=="R4",1); r1=find(retained=="R1",1);
    assert(~isempty(r4)&&~isempty(r1),'R4 or R1 was not retained in the full network.');
    distanceRows(end+1,:)=[expectedID(r4),expectedID(r1),2701];
end
writeRows([prefix '.dis'],distanceRows,'%d %d %.4f\n');
fprintf('Bundle distance observations: %d\n',size(distanceRows,1));
if numberOfBundleDistances==2
 fprintf('  F1--F6 = 2478 mm; R4--R1 = 2701 mm\n');
 fprintf('  Do NOT treat either as an independent check distance.\n');
else
 fprintf('  F1--F6 = 2478 mm (single scale distance).\n');
end
cfg=[1;0;0.001;0.0001;0.01;3.2;5;1;0.01];
writeRows([prefix '.cfg'],cfg,'%g\n');
writetable(table(retained,expectedID,'VariableNames',{'Target','ID'}),[prefix '_point_map.csv']);
fprintf('Prepared ONE nine-image adjustment: %s.*\n',prefix);
reply=input('Run the nine-image course adjustment? [y/N]: ','s');
if ~strcmpi(strtrim(reply),'y'),return;end
[eop,iop,con,tie,Cx,vr_img,vr_eop,vr_dis,vr_con]=full_matlab_bundle(prefix);
fprintf('Returned RMS: %.3f px (valid ONLY if course explicitly reports convergence)\n',sqrt(mean(vr_img(:,1).^2)));
save([prefix '_results.mat'],'distanceRows','numberOfBundleDistances','eop','iop','con','tie','Cx','vr_img','vr_eop','vr_dis','vr_con','E','TallNew','R','translation','summary');
end
function [X,rk]=triangulateWeightedLinear(obs,cameras,angles,IOP)
% Each observation gives two linearized collinearity equations:
% (m1 - xnorm*m3)*(X-C)=0; (m2 - ynorm*m3)*(X-C)=0
n=size(obs,1); AA=zeros(2*n,3); bb=zeros(2*n,1);
for k=1:n
    im=obs(k,2); M=rotW2C(angles(im,:)); C=cameras(im,:)';
    xn=-(obs(k,3)-IOP(1))/IOP(3);
    yn= (obs(k,4)-IOP(2))/IOP(3);
    G=[M(1,:)-xn*M(3,:);M(2,:)-yn*M(3,:)];
    AA(2*k-1:2*k,:)=G; bb(2*k-1:2*k)=G*C;
end
rk=rank(AA);X=(AA\bb)';
end
function errors=pointErrors(X,obs,cameras,angles,IOP)
errors=zeros(size(obs,1),1);
for k=1:size(obs,1)
  im=obs(k,2); proj=project(X,[cameras(im,:) angles(im,:)],IOP);
  errors(k)=norm(proj-obs(k,3:4));
end
end
function cost=pointCost(X,obs,cameras,angles,IOP)
e=pointErrors(X,obs,cameras,angles,IOP);
if any(~isfinite(e)) || max(abs(e))>1e7
 cost=1e18;
else
 cost=mean(e.^2);
end
end
function [R,t,rms]=rigidFit(initial,adjusted)
center0=mean(initial,1)';center1=mean(adjusted,1)';
H=(initial-center0')'*(adjusted-center1');
[U,~,V]=svd(H);
R=V*diag([1,1,det(V*U')])*U';
t=center1-R*center0;
d=(R*initial'+t)-adjusted';
rms=sqrt(mean(sum(d.^2,1)));
end
function ix=matchXYZ(src,dst,maxError,tag)
D=zeros(size(src,1),size(dst,1));
for j=1:size(src,1),D(j,:)=sqrt(sum((dst-src(j,:)).^2,2));end
[val,ix]=min(D,[],2);
assert(all(val<maxError),'Some %s tie coordinates have no matching all-nine initialization',tag);
assert(numel(unique(ix))==numel(ix),'Ambiguous %s tie matching',tag);
end
function M=rotW2C(wpk)
a=deg2rad(wpk(1)); b=deg2rad(wpk(2)); c=deg2rad(wpk(3));
Rx=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
Ry=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
Rz=[cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1];
M=(Rx*Ry*Rz)';
end
function ang=rot2WPK(M)
A=M';b=asin(max(-1,min(1,A(1,3))));
if abs(cos(b))>1e-7
 a=atan2(-A(2,3),A(3,3));c=atan2(-A(1,2),A(1,1));
else
 a=atan2(A(3,2),A(2,2)); c=0;
end
ang=rad2deg([a b c]);
end
function pix=project(X,p,IOP)
M=rotW2C(p(4:6));Q=(M*(X-p(1:3))')';
z=Q(:,3);z(abs(z)<1e-9)=NaN;
pix=[IOP(1)-IOP(3)*Q(:,1)./z,IOP(2)+IOP(3)*Q(:,2)./z];
end
function v=poseCost(x,X,obs,IOP)
p=[x(1:3)*1000 x(4:6)*30];pred=project(X,p,IOP);
if any(~isfinite(pred(:)))||any(abs(pred(:))>1e8),v=1e20;return;end
d=pred-obs;v=mean(sum(d.^2,2))+max(0,norm(p(1:3))-20000)^2/200;
end
function writeRows(file,A,fmt)
f=fopen(file,'w');assert(f>=0,'Cannot write output');c=onCleanup(@()fclose(f));
for k=1:size(A,1),fprintf(f,fmt,A(k,:));end
end
