ods rtf style=styles.Myway1
file="E:OUTPUT.rtf";

*Advanced Regression Analysis

ods html style=styles.Myway1;

*a) Ordinal Logistic Regression;
proc import out=stroke
	file="E:\SASdata\stroke.sav"
	dbms=spss replace;
run;
proc sort data=stroke out=stroke;
	by rankin0 anticlot smoker af gender age;
run;
ods graphics on / ANTIALIASMAX=12200 noborder;
proc logistic data=stroke order=data plots=all;
	title Ordinal Logistic Regression;
	class rankin0 gender(ref='Female') af(ref='No')
		anticlot(ref='None') smoker(ref='No') / param=ref;
	model rankin0 = age gender af smoker anticlot / rsq stb expb clparm=wald
		link=logit  aggregate scale=none;
	output out=c_probs p=c_probs predprobs=individual;
	oddsratio age / diff=ref;
	oddsratio gender / diff=ref;
	oddsratio af / diff=ref;
	oddsratio smoker / diff=ref;
	oddsratio anticlot / diff=ref;
	effectplot interaction(plotby=rankin0 x=gender) / clm;
	effectplot interaction(plotby=rankin0 x=af) / clm;
	effectplot interaction(plotby=rankin0 x=anticlot) / clm;
	effectplot interaction(plotby=rankin0 x=smoker) / clm;
	effectplot interaction(plotby=rankin0 x=gender) / individual;
	effectplot interaction(plotby=rankin0 x=af) / individual;
	effectplot interaction(plotby=rankin0 x=anticlot) / individual;
	effectplot interaction(plotby=rankin0 x=smoker) / individual;
	effectplot interaction(sliceby=rankin0 x=smoker) / individual polybar;
	effectplot interaction(sliceby=rankin0 x=af) / individual polybar;
	effectplot interaction(sliceby=rankin0 x=anticlot) / individual polybar;
	effectplot interaction(sliceby=rankin0 x=gender) / individual polybar;
run;

*Regressor 1: Age;
proc sgplot data=c_probs;
	title 'Regressor 1: Age';
	loess x=age y=IP_no_symptoms
	/ markerattrs=(symbol=circle) lineattrs=(thickness=3);
	yaxis label='Cumulative P of No Symptoms';
run;

*Regressor 2: Gender;
data female; set C_probs(where=(gender=1));
	IP_no_symptoms_female = IP_no_symptoms;
run;
data male; set C_probs(where=(gender=0));
	IP_no_symptoms_male = IP_no_symptoms;
run;
data new; set female male; run;
proc sgplot data=new;
	title 'Regressor 2: Gender';
	histogram IP_no_symptoms_female
		/ legendlabel="Females" name='a';
	histogram IP_no_symptoms_male
		/ legendlabel="Males" name='b';
	density IP_no_symptoms_female
		/ lineattrs=(color='black' thickness=3);
	density IP_no_symptoms_male 
		/ lineattrs=(color='black' thickness=3);
	keylegend 'a' 'b';
	xaxis label='Cumulative P of No Symptoms';
run;

*Regressor 3: Atrial Fibrillation;
data af1; set C_probs(where=(af=1));
	IP_no_symptoms_af1 = IP_no_symptoms;
run;
data af0; set C_probs(where=(af=0));
	IP_no_symptoms_af0 = IP_no_symptoms;
run;
data new; set af1 af0; run;
proc sgplot data=new;
	title 'Regressor 3: Atrial Fibrillation';
	histogram IP_no_symptoms_af1
		/ legendlabel="Atrial Fibrillation" name='b';
	histogram IP_no_symptoms_af0
		/ legendlabel="Healthy Pulse" name='a';
	density IP_no_symptoms_af0
		/ lineattrs=(color='black' thickness=3);
	density IP_no_symptoms_af1
		/ lineattrs=(color='black' thickness=3);
	keylegend 'a' 'b';
	xaxis label='Cumulative P of No Symptoms';
run;

*Regressor 4: Smokers;
data smoker; set C_probs(where=(smoker=1));
	IP_no_symptoms_smoke = IP_no_symptoms;
run;
data nonsmoker; set C_probs(where=(smoker=0));
	IP_no_symptoms_nonsmoke = IP_no_symptoms;
run;
data new; set smoker nonsmoker; run;
proc sgplot data=new;
	title 'Regressor 4: Smokers';
	histogram IP_no_symptoms_smoke
		/ legendlabel="Smokers" name='a';
	histogram IP_no_symptoms_nonsmoke
		/ legendlabel="Non-Smokers" name='b';
	density IP_no_symptoms_smoke
		/ lineattrs=(color='black' thickness=3);
	density IP_no_symptoms_nonsmoke
		/ lineattrs=(color='black' thickness=3);
	keylegend 'a' 'b';
	xaxis label='Cumulative P of No Symptoms';
run;

*Regressor 5: Anti-Clot Medication;
data ACM0; set C_probs(where=(anticlot=0));
	IP_no_symptoms_ACM0 = IP_no_symptoms;
run;
data ACM1; set C_probs(where=(anticlot=1));
	IP_no_symptoms_ACM1 = IP_no_symptoms;
run;
data ACM2; set C_probs(where=(anticlot=2));
	IP_no_symptoms_ACM2 = IP_no_symptoms;
run;
data ACM3; set C_probs(where=(anticlot=3));
	IP_no_symptoms_ACM3 = IP_no_symptoms;
run;
data new; set ACM0 ACM1 ACM2 ACM3; run;

proc sgplot data=new;
	title 'Regressor 5: Anti-Clot Medication';
	histogram IP_no_symptoms_ACM3
		/ legendlabel="Warfarin" name='d';
	density IP_no_symptoms_ACM3
		/ lineattrs=(color='black' thickness=3);
	histogram IP_no_symptoms_ACM2
		/ legendlabel="Heparin" name='c';
	density IP_no_symptoms_ACM2
		/ lineattrs=(color='black' thickness=3);
	histogram IP_no_symptoms_ACM1
		/ legendlabel="Aspirin" name='b';
	density IP_no_symptoms_ACM1
		/ lineattrs=(color='black' thickness=3);
	histogram IP_no_symptoms_ACM0
		/ legendlabel="None" name='a';
	density IP_no_symptoms_ACM0
		/ lineattrs=(color='black' thickness=3);
	keylegend 'a' 'b' 'c' 'd';
	xaxis label='Cumulative P of No Symptoms';
run;

proc sgplot data=new;
	title 'Regressor 5: Anti-Clot Medication';
	density IP_no_symptoms_ACM3	/ lineattrs=(thickness=3)
	 	legendlabel="Warfarin" name='d' SCALE=percent;
	density IP_no_symptoms_ACM2 / lineattrs=(thickness=3)
		legendlabel="Heparin" name='c' SCALE=percent;
	density IP_no_symptoms_ACM1	/ lineattrs=(thickness=3)
		legendlabel="Aspirin" name='b' SCALE=percent;
	density IP_no_symptoms_ACM0	/ lineattrs=(thickness=3)
		legendlabel="None" name='a' SCALE=percent;
	keylegend 'a' 'b' 'c' 'd';
	xaxis label='Cumulative P of No Symptoms';
run;

*lots of warfarin users have af?;
proc freq data=stroke;
	table anticlot*af/ norow nopercent ;
run;
proc freq data=stroke;
	table rankin0*af/ norow nopercent ;
run;
proc sgpanel data=c_probs;
	panelby af;
	hbar  anticlot;
	hline anticlot;
run;
title ; ods graphics off;


*b)Multinomial Logistic Regression;
proc import out=stroke
	file="E:\SASdata\stroke.sav"
	dbms=spss replace;
run;
proc sort data=stroke out=stroke;
	by rankin0 anticlot smoker af gender age;
run;
proc format library=work.formats;
	value newrank	0='No Disability'
					1='Moderate Dis.'
					2='Severe Disab.';
run;
data stroke_b; set stroke;
	if rankin0=0 or rankin0=1 then rankin1=0;
	if rankin0=2 or rankin0=3 then rankin1=1;
	if rankin0=4 or rankin0=5 then rankin1=2;
	format rankin1 newrank.;
run;
proc freq data=stroke_b;
	table rankin0*rankin1 / nocol norow nopercent;
run;
ods graphics on / ANTIALIASMAX=12200 noborder;
proc logistic data=stroke_b order=data plots=all OUTEST=ests;
	title Ordinal Logistic Regression;
	class rankin1 gender(ref='Female') af(ref='No')
		anticlot(ref='None') smoker(ref='No') / param=ref;
	model rankin1 = age gender af smoker anticlot / rsq stb expb clparm=wald
		link=glogit aggregate scale=none;
	output out=c_probs p=c_probs predprobs=individual;
	oddsratio age / diff=ref;
	oddsratio gender / diff=ref;
	oddsratio af / diff=ref;
	oddsratio smoker / diff=ref;
	oddsratio anticlot / diff=ref;
	effectplot interaction(plotby=rankin1 x=gender) / clm;
	effectplot interaction(plotby=rankin1 x=af) / clm;
	effectplot interaction(plotby=rankin1 x=anticlot) / clm;
	effectplot interaction(plotby=rankin1 x=smoker) / clm;
run;
title ; ods graphics off;


* Residuals
proc import out=residuals
	file="E:\resfil2b.sav"
	dbms=spss replace;
run;
ods graphics on / noborder height=10cm;
proc sgplot data=residuals;
	title Level-2 Residuals;
	scatter x=chipct y=mdist;
	xaxis label='Chi-Square Quantiles';
	yaxis label='Mahalanobis Distance';
run;
title ; ods graphics off;


* Survival Analysis
proc import out=stroke_survivalb
	file="E:\SASdata\stroke_survivalb.sav"
	dbms=spss replace;
run;

*a) Summary Statistics;
proc format library=work.formats;
	value alive		0='Survived'
					1='Dropped';
run;
data stroke_survivalb2;
	set stroke_survivalb;
	if time1<1500 then died = 1;
		else died=0;
	if time1<1500 then deadhist = age;
	else alivehist = age;
	format died alive.;
run;
proc sort data=stroke_survivalb2 out=stroke_survivalb2;
	by died;
run;
proc means data=stroke_survivalb2 maxdec=1;
	by died;
	var time1 age;
run;

ods graphics on / noborder height=8cm;
proc sgplot data=stroke_survivalb2(where=(died=1));
	title 'Variable Histogram: Survival Time';
	histogram time1;
	density time1 / type=kernel;
	xaxis label='Survival Time';
run;
proc freq data=stroke_survivalb2;
	table rehab * died / nopercent norow;
	table surgery * died / nopercent norow;
run;
title ; ods graphics off;

*b) Kaplan-Meier Rankin0;
ods graphics on / noborder height=13cm;
proc lifetest data=stroke_survivalb plots=(s);
	time time1*event(0);
	strata rankin0;
run;
ods graphics off;

*c) Kaplan-Meier Rehab;
ods graphics on / noborder height=13cm;
proc lifetest data=stroke_survivalb plots=(s);
	time time1*event(0);
	strata rehab;
run;
ods graphics off;

*d) Cox regression Rankin0 and Rehab;
 ods graphics on / noborder height=13cm;
proc phreg data=stroke_survivalb plots=(cumhaz survival ) COVS ;
	class rankin0 rehab surgery;
	model time1*event(0) = rankin0 rehab surgery;
	assess ph / RESAMPLE crpanel ;
	output out=temp resmart=Mresids resdev=Dresids ressch=Sresids;
run;
ods graphics off;

 ods graphics on / noborder height=13cm;
proc phreg data=relapse plots=(cumhaz survival) COVS(AGGREGATE);
	class chemo;
	model time*status(0)= chemo;
output out=temp resmart=Mresids
resdev=Dresids ressch=Sresids;
run;
ods graphics off;




ods rtf close;
