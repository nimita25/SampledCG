Algo = [3.0128,
3,
3.01,
3.589,
3];

SDPNAL = [12.91;
11.63;
12.98;
13.5;
13.17;
11.5;
14.35;
12.86;
12.71;
77.25;
65.22;
65.5;
79.91;
65.17;
83.37;
88.36;
82.34;
85.53;
20.11;
20.42;
20.24;
19.85;
20.3;
155.23;
171.66;
177.57;
21.65;
21.4;
21.55;
20.73;
456.57;
625.07;
747.75];


SDPT3 = [9.786;
9.79;
9.79;
9.79;
9.79;
9.78;
9.78;
9.78;
9.78;
61.07;
61.07;
61.05;
61.06;
61.11;
61.06;
61.05;
61.05;
61.05;
15.27;
15.27;
15.27;
15.27;
15.27;
137.34;
137.36;
137.36;
15.28;
15.27;
15.27;
15.27;
565.14;
381.53;
1107.8;
747.79];

SeDuMi = [9.812,
77.2,
15.32,
137.424,
381.69];

CGAL = [1;
2.11;
0.98;
0.98;
0.99;
1.23;
2.05;
0.72;
1.25;
1.31;
1.45;
1.11;
1.1;
1.13;
0.92;
1.26;
0.99;
0.95;
1.01;
0.86;
0.9;
1.1;
0.64;
0.94;
0.89;
0.94;
0.88;
0.69;
3.88;
0.69];

Vertices = [800;
800;
800;
800;
800;
800;
800;
800;
800;
2000;
2000;
2000;
2000;
2000;
2000;
2000;
2000;
2000;
1000;
1000;
1000;
1000;
1000;
3000;
3000;
3000;
1000;
1000;
1000;
1000;
5000;
5000;
7000;
7000];

p = figure;
semilogy(Vertices(1:length(Algo)), Algo,'o','MarkerSize',10)
title('Comparison of Memory Requirement (in MB)');
hold on;
semilogy(Vertices(1:length(SDPNAL)), SDPNAL,'+','MarkerSize',10)
semilogy(Vertices(1:length(SDPT3)), SDPT3,'s','MarkerSize',10)
semilogy(Vertices(1:length(SeDuMi)), SeDuMi,'d','MarkerSize',10)
semilogy(Vertices(1:length(CGAL)), CGAL,'^','MarkerSize',10)
hold off;
legend({'Algorithm 2','SDPNAL+','SDPT3','SeDuMi','CGAL (R=10)'},'Location','southeast');
ylabel('Memory required in MB');
xlabel('Number of vertices');
%saveas(p,'MemoryPlot.png');