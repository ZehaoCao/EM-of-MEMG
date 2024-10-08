%Generate forecast output of DERs: scripted by Zehao Cao
function [windmax,windmin,pvmax,pvmin] = DERs_org_data()
DER_factor = 1;
wind_pow = 0.14;
wind_pow = wind_pow*1.5*DER_factor;
windmax = wind_pow*[ 																																															
  101.1	102.4 101.4	93.2	98.6 103.6	109.2	104	101.4	61.6	65.2	70.4	73	66.5	67.8	72.8	106.6	98.3	104.2	91.26	98.8	106.6	109.2	112.3
    ];																								
windmin = wind_pow*[																																															
130	130	101.4	101.4	130	119.6 130	130	101.4	83.2	130	119.6 130	130	101.4	83.2	130	119.6 130	130	101.4	130	130	119.6
    ];																																																
pv_pow = 0.2;
pv_pow = pv_pow*1.5*DER_factor;
pvmax = pv_pow*[ 																																															
 0 	0 	0 	0 	0 	4.068 	5.566 	6.211 	7.773 	7.714 	8.022 	8.746 	7.215 	7.038 	6.704 	6.445 	5.397 	2.877 	1.171 	0 	0 	0 	0 	0 
    ];																								
pvmin = 0.9*pv_pow*[																																														
 0 	0 	0 	0 	0 	1.068 	2.566 	3.211 	4.773 	4.714 	5.022 	4.746 	4.215 	5.338 	3.704 	4.445 	3.497 	1.877 	1.171 	0 	0 	0 	0 	0
    ];
end																								
