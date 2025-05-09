reset
set style data lines
set nokey

set xlabel 'x (m)'
set title 'Alpha\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:2,\
 "./datasets/result_CPU1_TIME0.out" u 1:2,\
 "./datasets/result_CPU2_TIME0.out" u 1:2,\
 "./datasets/result_CPU0_TIME1.out" u 1:2,\
 "./datasets/result_CPU1_TIME1.out" u 1:2,\
 "./datasets/result_CPU2_TIME1.out" u 1:2
pause(-1)

set title 'Density\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:3,\
 "./datasets/result_CPU1_TIME0.out" u 1:3,\
 "./datasets/result_CPU2_TIME0.out" u 1:3,\
 "./datasets/result_CPU0_TIME1.out" u 1:3,\
 "./datasets/result_CPU1_TIME1.out" u 1:3,\
 "./datasets/result_CPU2_TIME1.out" u 1:3
pause(-1)

set title 'Pressure\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:4,\
 "./datasets/result_CPU1_TIME0.out" u 1:4,\
 "./datasets/result_CPU2_TIME0.out" u 1:4,\
 "./datasets/result_CPU0_TIME1.out" u 1:4,\
 "./datasets/result_CPU1_TIME1.out" u 1:4,\
 "./datasets/result_CPU2_TIME1.out" u 1:4
pause(-1)

set title 'Cobase\_11\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:5,\
 "./datasets/result_CPU1_TIME0.out" u 1:5,\
 "./datasets/result_CPU2_TIME0.out" u 1:5,\
 "./datasets/result_CPU0_TIME1.out" u 1:5,\
 "./datasets/result_CPU1_TIME1.out" u 1:5,\
 "./datasets/result_CPU2_TIME1.out" u 1:5
pause(-1)

set title 'Cobase\_12\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:6,\
 "./datasets/result_CPU1_TIME0.out" u 1:6,\
 "./datasets/result_CPU2_TIME0.out" u 1:6,\
 "./datasets/result_CPU0_TIME1.out" u 1:6,\
 "./datasets/result_CPU1_TIME1.out" u 1:6,\
 "./datasets/result_CPU2_TIME1.out" u 1:6
pause(-1)

set title 'Cobase\_13\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:7,\
 "./datasets/result_CPU1_TIME0.out" u 1:7,\
 "./datasets/result_CPU2_TIME0.out" u 1:7,\
 "./datasets/result_CPU0_TIME1.out" u 1:7,\
 "./datasets/result_CPU1_TIME1.out" u 1:7,\
 "./datasets/result_CPU2_TIME1.out" u 1:7
pause(-1)

set title 'Cobase\_21\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:8,\
 "./datasets/result_CPU1_TIME0.out" u 1:8,\
 "./datasets/result_CPU2_TIME0.out" u 1:8,\
 "./datasets/result_CPU0_TIME1.out" u 1:8,\
 "./datasets/result_CPU1_TIME1.out" u 1:8,\
 "./datasets/result_CPU2_TIME1.out" u 1:8
pause(-1)

set title 'Cobase\_22\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:9,\
 "./datasets/result_CPU1_TIME0.out" u 1:9,\
 "./datasets/result_CPU2_TIME0.out" u 1:9,\
 "./datasets/result_CPU0_TIME1.out" u 1:9,\
 "./datasets/result_CPU1_TIME1.out" u 1:9,\
 "./datasets/result_CPU2_TIME1.out" u 1:9
pause(-1)

set title 'Cobase\_23\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:10,\
 "./datasets/result_CPU1_TIME0.out" u 1:10,\
 "./datasets/result_CPU2_TIME0.out" u 1:10,\
 "./datasets/result_CPU0_TIME1.out" u 1:10,\
 "./datasets/result_CPU1_TIME1.out" u 1:10,\
 "./datasets/result_CPU2_TIME1.out" u 1:10
pause(-1)

set title 'Cobase\_31\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:11,\
 "./datasets/result_CPU1_TIME0.out" u 1:11,\
 "./datasets/result_CPU2_TIME0.out" u 1:11,\
 "./datasets/result_CPU0_TIME1.out" u 1:11,\
 "./datasets/result_CPU1_TIME1.out" u 1:11,\
 "./datasets/result_CPU2_TIME1.out" u 1:11
pause(-1)

set title 'Cobase\_32\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:12,\
 "./datasets/result_CPU1_TIME0.out" u 1:12,\
 "./datasets/result_CPU2_TIME0.out" u 1:12,\
 "./datasets/result_CPU0_TIME1.out" u 1:12,\
 "./datasets/result_CPU1_TIME1.out" u 1:12,\
 "./datasets/result_CPU2_TIME1.out" u 1:12
pause(-1)

set title 'Cobase\_33\_SG\_aluminum'
plot "./datasets/result_CPU0_TIME0.out" u 1:13,\
 "./datasets/result_CPU1_TIME0.out" u 1:13,\
 "./datasets/result_CPU2_TIME0.out" u 1:13,\
 "./datasets/result_CPU0_TIME1.out" u 1:13,\
 "./datasets/result_CPU1_TIME1.out" u 1:13,\
 "./datasets/result_CPU2_TIME1.out" u 1:13
pause(-1)

set title 'Density\_Mixture'
plot "./datasets/result_CPU0_TIME0.out" u 1:14,\
 "./datasets/result_CPU1_TIME0.out" u 1:14,\
 "./datasets/result_CPU2_TIME0.out" u 1:14,\
 "./datasets/result_CPU0_TIME1.out" u 1:14,\
 "./datasets/result_CPU1_TIME1.out" u 1:14,\
 "./datasets/result_CPU2_TIME1.out" u 1:14
pause(-1)

set title 'Pressure\_Mixture'
plot "./datasets/result_CPU0_TIME0.out" u 1:15,\
 "./datasets/result_CPU1_TIME0.out" u 1:15,\
 "./datasets/result_CPU2_TIME0.out" u 1:15,\
 "./datasets/result_CPU0_TIME1.out" u 1:15,\
 "./datasets/result_CPU1_TIME1.out" u 1:15,\
 "./datasets/result_CPU2_TIME1.out" u 1:15
pause(-1)

set title 'Total\_energy\_Mixture'
plot "./datasets/result_CPU0_TIME0.out" u 1:16,\
 "./datasets/result_CPU1_TIME0.out" u 1:16,\
 "./datasets/result_CPU2_TIME0.out" u 1:16,\
 "./datasets/result_CPU0_TIME1.out" u 1:16,\
 "./datasets/result_CPU1_TIME1.out" u 1:16,\
 "./datasets/result_CPU2_TIME1.out" u 1:16
pause(-1)

set title 'Stress\_tensor\_11'
plot "./datasets/result_CPU0_TIME0.out" u 1:17,\
 "./datasets/result_CPU1_TIME0.out" u 1:17,\
 "./datasets/result_CPU2_TIME0.out" u 1:17,\
 "./datasets/result_CPU0_TIME1.out" u 1:17,\
 "./datasets/result_CPU1_TIME1.out" u 1:17,\
 "./datasets/result_CPU2_TIME1.out" u 1:17
pause(-1)

set title 'Stress\_tensor\_12'
plot "./datasets/result_CPU0_TIME0.out" u 1:18,\
 "./datasets/result_CPU1_TIME0.out" u 1:18,\
 "./datasets/result_CPU2_TIME0.out" u 1:18,\
 "./datasets/result_CPU0_TIME1.out" u 1:18,\
 "./datasets/result_CPU1_TIME1.out" u 1:18,\
 "./datasets/result_CPU2_TIME1.out" u 1:18
pause(-1)

set title 'Velocity\_x'
plot "./datasets/result_CPU0_TIME0.out" u 1:19,\
 "./datasets/result_CPU1_TIME0.out" u 1:19,\
 "./datasets/result_CPU2_TIME0.out" u 1:19,\
 "./datasets/result_CPU0_TIME1.out" u 1:19,\
 "./datasets/result_CPU1_TIME1.out" u 1:19,\
 "./datasets/result_CPU2_TIME1.out" u 1:19
pause(-1)

set title 'Velocity\_y'
plot "./datasets/result_CPU0_TIME0.out" u 1:20,\
 "./datasets/result_CPU1_TIME0.out" u 1:20,\
 "./datasets/result_CPU2_TIME0.out" u 1:20,\
 "./datasets/result_CPU0_TIME1.out" u 1:20,\
 "./datasets/result_CPU1_TIME1.out" u 1:20,\
 "./datasets/result_CPU2_TIME1.out" u 1:20
pause(-1)

set title 'Velocity\_Mixture'
plot "./datasets/result_CPU0_TIME0.out" u 1:21,\
 "./datasets/result_CPU1_TIME0.out" u 1:21,\
 "./datasets/result_CPU2_TIME0.out" u 1:21,\
 "./datasets/result_CPU0_TIME1.out" u 1:21,\
 "./datasets/result_CPU1_TIME1.out" u 1:21,\
 "./datasets/result_CPU2_TIME1.out" u 1:21
pause(-1)

set title 'AMR level'
plot "./datasets/result_CPU0_TIME0.out" u 1:22,\
 "./datasets/result_CPU1_TIME0.out" u 1:22,\
 "./datasets/result_CPU2_TIME0.out" u 1:22,\
 "./datasets/result_CPU0_TIME1.out" u 1:22,\
 "./datasets/result_CPU1_TIME1.out" u 1:22,\
 "./datasets/result_CPU2_TIME1.out" u 1:22
pause(-1)

set title 'Xi'
plot "./datasets/result_CPU0_TIME0.out" u 1:23,\
 "./datasets/result_CPU1_TIME0.out" u 1:23,\
 "./datasets/result_CPU2_TIME0.out" u 1:23,\
 "./datasets/result_CPU0_TIME1.out" u 1:23,\
 "./datasets/result_CPU1_TIME1.out" u 1:23,\
 "./datasets/result_CPU2_TIME1.out" u 1:23
pause(-1)

