%initialization CobraToolbox and GUROBI as solver
initCobraToolbox
changeCobraSolver('gurobi','all')
setRavenSolver('gurobi')
%% Load the model

getFullPath(fullfile(cd, 'Bacteroides_thetaiotaomicron_VPI_5482.xml'))
cd [your_path];
model=readCbModel('Bacteroides_thetaiotaomicron_VPI_5482','fileType','SBML');

%% Rename the model
model_BT = model;

%% Curation 
%Update metabolites formula
pos1 = strcmp(model_BT.metNames,'L-Cystathionine');
k= find(pos1);
model_BT.metFormulas(k);
model_BT.metFormulas(k)= {'C7H14N2O4S'};

pos2 = strcmp(model_BT.metNames,'Pyridoxal 5''-phosphate');
k2= find(pos2);
model_BT.metFormulas(k2);
model_BT.metFormulas(k2)= {'C8H8NO6P'};

pos3 = strcmp(model_BT.metNames,'Oxidized ferredoxin');
k3= find(pos3);
model_BT.metFormulas(k3);
model_BT.metFormulas(k3)= {'Fe2S2X'};

pos4 = strcmp(model_BT.metNames,'2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate');
k4= find(pos4);
model_BT.metFormulas(k4);
model_BT.metFormulas(k4)= {'C9H12N5O13P3'};

% Update existing reactions
model_BT = removeRxns(model_BT, {'R_FRDO'});
model_BT=addReaction(model_BT, 'R_FRDO','reactionName', 'FRDO', 'reactionformula', 'nadp[C_c] + fdxrd[C_c] + h[C_c] <=> nadph[C_c] + fdxox[C_c]','geneRule','( NP_811748_1 and NP_811749_1 )','lowerBound',-1000,'upperBound',1000);%original adp[C_c] + fdxrd[C_c] 	<=>	h[C_c] + nadph[C_c] + fdxox[C_c]

model_BT = removeRxns(model_BT, {'R_PANTS'});
model_BT=addReaction(model_BT, 'R_PANTS','reactionName', 'Pantothenate synthase', 'reactionformula', 'atp[C_c] + ala_B[C_c] + pant__R[C_c] -> 2 h[C_c] + amp[C_c] + ppi[C_c] + pnto__R[C_c]','geneRule','NP_813219_1','lowerBound',0,'upperBound',1000);%original atp[C_c] + ala_B[C_c] + pant__R[C_c] 	->	h[C_c] + amp[C_c] + ppi[C_c] + pnto__R[C_c]

model_BT = removeRxns(model_BT, {'R_PTPATi'});
model_BT=addReaction(model_BT, 'R_PTPATi','reactionName', 'Pantetheine-phosphate adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + pan4p[C_c] 	->	ppi[C_c] + dpcoa[C_c]','geneRule','NP_811946_1','lowerBound',0,'upperBound',1000);%original atp[C_c] + h[C_c] + pan4p[C_c] 	->	ppi[C_c] + dpcoa[C_c]

model_BT = removeRxns(model_BT, {'R_PUNP4'});
model_BT=addReaction(model_BT, 'R_PUNP4','reactionName', 'Purine-nucleoside phosphorylase (Deoxyguanosine)', 'reactionformula', 'pi[C_c] + dgsn[C_c] 	<=>	2dr1p[C_c] + gua[C_c] + h[C_c]','geneRule','( NP_810794_1 or NP_813465_1 )','lowerBound',-1000,'upperBound',1000);%original pi[C_c] + dgsn[C_c] 	<=>	2dr1p[C_c] + gua[C_c]

model_BT = removeRxns(model_BT, {'R_PNP'});
model_BT=addReaction(model_BT, 'R_PNP','reactionName', 'Purine nucleoside phosphorylase', 'reactionformula', 'pi[C_c] + rnam[C_c] <=>	ncam[C_c] + r1p[C_c]','geneRule','( NP_810794_1 or NP_813465_1 )','lowerBound',-1000,'upperBound',1000);%original pi[C_c] + rnam[C_c] 	<=>	h[C_c] + ncam[C_c] + r1p[C_c]

model_BT = removeRxns(model_BT, {'R_MTHFR3'});
model_BT=addReaction(model_BT, 'R_MTHFR3','reactionName', '5 10 methylenetetrahydrofolatereductase  NADPH', 'reactionformula', 'h[C_c] + nadph[C_c] + mlthf[C_c] ->	nadp[C_c] + 5mthf[C_c]','geneRule','NP_812732_1','lowerBound',0,'upperBound',1000);%original 2 h[C_c] + nadph[C_c] + mlthf[C_c] 	->	nadp[C_c] + 5mthf[C_c]

model_BT = removeRxns(model_BT, {'R_NADDP'});
model_BT=addReaction(model_BT, 'R_NADDP','reactionName', 'NAD diphosphatase', 'reactionformula', 'h2o[C_c] + nad[C_c] 	->	3 h[C_c] + amp[C_c] + nmn[C_c]','geneRule','( NP_809807_1 or NP_810457_1 )','lowerBound',0,'upperBound',1000);%original h2o[C_c] + nad[C_c] 	->	2 h[C_c] + amp[C_c] + nmn[C_c]

model_BT = removeRxns(model_BT, {'R_OACT'});
model_BT=addReaction(model_BT, 'R_OACT','reactionName', 'O-antigen Acetyl-Transferase', 'reactionformula', '2 h[C_c] + accoa[C_c] + udcpo4[C_c] -> coa[C_c] + udcpo5[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original accoa[C_c] + udcpo4[C_c] 	->	coa[C_c] + udcpo5[C_c]

model_BT = removeRxns(model_BT, {'R_ALLTN'});
model_BT=addReaction(model_BT, 'R_ALLTN','reactionName', 'Allantoinase', 'reactionformula', 'h2o[C_c] + alltn[C_c] -> alltt[C_c]','geneRule','NP_809163_1','lowerBound',0,'upperBound',1000);%original h2o[C_c] + alltn[C_c] 	->	h[C_c] + alltt[C_c]

model_BT = removeRxns(model_BT, {'R_IG3PS_1'});
model_BT=addReaction(model_BT, 'R_IG3PS_1','reactionName', 'Allantoinase', 'reactionformula', 'gln__L[C_c] + prlp[C_c] -> glu__L[C_c] + h[C_c] + aicar[C_c] + eig3p[C_c]','geneRule','NP_810291_1','lowerBound',0,'upperBound',1000);%original gln__L[C_c] + prlp[C_c] 	->	glu__L[C_c] + 2 h[C_c] + aicar[C_c] + eig3p[C_c]

model_BT = removeRxns(model_BT, {'R_FMNRx'});
model_BT=addReaction(model_BT, 'R_FMNRx','reactionName', 'FMN reductase', 'reactionformula', '2 h[C_c] + nadh[C_c] + fmn[C_c] -> nad[C_c] + fmnh2[C_c]','geneRule','NP_810980_1','lowerBound',0,'upperBound',1000);%original h[C_c] + nadh[C_c] + fmn[C_c] 	->	nad[C_c] + fmnh2[C_c]

model_BT = removeRxns(model_BT, {'R_FMNRy'});
model_BT=addReaction(model_BT, 'R_FMNRy','reactionName', 'FMN reductase  NADPH dependent', 'reactionformula', '2 h[C_c] + nadph[C_c] + fmn[C_c] -> nadp[C_c] + fmnRD[C_c]','geneRule','NP_811057_1','lowerBound',0,'upperBound',1000);%original h[C_c] + nadph[C_c] + fmn[C_c] 	->	nadp[C_c] + fmnRD[C_c]

model_BT = removeRxns(model_BT, {'R_PRAGSr'});
model_BT=addReaction(model_BT, 'R_PRAGSr','reactionName', 'Phosphoribosylglycinamide synthase', 'reactionformula', 'atp[C_c] + gly[C_c] + pram[C_c] <=>	adp[C_c] + pi[C_c] + 2 h[C_c] + gar[C_c]','geneRule','NP_811057_1','lowerBound',-1000,'upperBound',1000);%original atp[C_c] + gly[C_c] + pram[C_c] 	<=>	adp[C_c] + pi[C_c] + h[C_c] + gar[C_c]

model_BT = removeRxns(model_BT, {'R_NMNAT'});
model_BT=addReaction(model_BT, 'R_NMNAT','reactionName', 'Nicotinamide-nucleotide adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + nmn[C_c] ->	nad[C_c] + ppi[C_c]','geneRule','( NP_810924_1 or NP_811727_1 )','lowerBound',0,'upperBound',1000);%original atp[C_c] + h[C_c] + nmn[C_c] 	->	nad[C_c] + ppi[C_c]

model_BT = removeRxns(model_BT, {'R_GUAD'});
model_BT=addReaction(model_BT, 'R_GUAD','reactionName', 'Guanine deaminase', 'reactionformula', 'h2o[C_c] + h[C_c] + gua[C_c] -> nh3[C_c] + xan[C_c]','geneRule','NP_809234_1','lowerBound',0,'upperBound',1000);%original fmnh2[C_c] + 2 salchs4fe[C_c] 	->	2 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 salchs4[C_c]

model_BT = removeRxns(model_BT, {'R_FEENTERR2'});
model_BT=addReaction(model_BT, 'R_FEENTERR2','reactionName', 'Fe-enterobactin reduction (Fe(III)-unloading)', 'reactionformula', 'fmnh2[C_c] + 2 feenter[C_c] -> 3 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 enter[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original fmnh2[C_c] + 2 feenter[C_c] 	->	2 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 enter[C_c]

model_BT = removeRxns(model_BT, {'R_GUAPRT'});
model_BT=addReaction(model_BT, 'R_GUAPRT','reactionName', 'Guanine phosphoribosyltransferase', 'reactionformula', 'h[C_c] + prpp[C_c] + gua[C_c] ->	ppi[C_c] + gmp[C_c]','geneRule','NP_813297_1','lowerBound',0,'upperBound',1000);%original prpp[C_c] + gua[C_c] 	->	ppi[C_c] + gmp[C_c]

model_BT = removeRxns(model_BT, {'R_FMNAT'});
model_BT=addReaction(model_BT, 'R_FMNAT','reactionName', 'FMN adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + fmn[C_c] -> ppi[C_c] + fad[C_c]','geneRule','NP_811456_1','lowerBound',0,'upperBound',1000);%original atp[C_c] + h[C_c] + fmn[C_c] 	->	ppi[C_c] + fad[C_c]

model_BT = removeRxns(model_BT, {'R_3NUCLE4'});
model_BT=addReaction(model_BT, 'R_3NUCLE4','reactionName', '3  nucleotidase  guanosine 3  phosphate', 'reactionformula', 'h2o[C_c] + 3gmp[C_c] 	-> h[C_c] + pi[C_c] + gsn[C_c]','geneRule','NP_810149_1','lowerBound',0,'upperBound',1000);%original h2o[C_c] + 3gmp[C_c] 	->	pi[C_c] + gsn[C_c]

model_BT = removeRxns(model_BT, {'R_NTD9'});
model_BT=addReaction(model_BT, 'R_NTD9','reactionName', '5-nucleotidase (GMP)', 'reactionformula', 'h2o[C_c] + gmp[C_c] ->	h[C_c] + pi[C_c] + gsn[C_c]','geneRule','NP_810149_1','lowerBound',0,'upperBound',1000);%original h2o[C_c] + gmp[C_c] 	->	pi[C_c] + gsn[C_c]

model_BT = removeRxns(model_BT, {'R_SO3R'});
model_BT=addReaction(model_BT, 'R_SO3R','reactionName', 'Sulfite reductase', 'reactionformula', '4 h[C_c] + 3 nadh[C_c] + so3[C_c]  -> 3 h2o[C_c] + 3 nad[C_c] + h2s[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original 5 h[C_c] + 3 nadh[C_c] + so3[C_c]  -> 3 h2o[C_c] + 3 nad[C_c] + h2s[C_c] 

model_BT = removeRxns(model_BT, {'R_PAPSR'});
model_BT=addReaction(model_BT, 'R_PAPSR','reactionName', 'Phosphoadenylyl-sulfate reductase (thioredoxin)', 'reactionformula', 'paps[C_c] + trdrd[C_c]  -> h[C_c] + pap[C_c] + so3[C_c] + trdox[C_c]','geneRule','( NP_809131_1 or NP_811142_1 or NP_811224_1 )','lowerBound',0,'upperBound',1000);%original paps[C_c] + trdrd[C_c]  -> 2 h[C_c] + pap[C_c] + so3[C_c] + trdox[C_c] 

model_BT = removeRxns(model_BT, {'R_URCN'});
model_BT=addReaction(model_BT, 'R_URCN','reactionName', 'Urocanase', 'reactionformula', 'h2o[C_c] + urcan[C_c] -> 4izp[C_c] + h[C_c]','geneRule','NP_811606_1','lowerBound',0,'upperBound',1000);%original h2o[C_c] + urcan[C_c] 	->	4izp[C_c]

model_BT = removeRxns(model_BT, {'R_RBFSb'});
model_BT=addReaction(model_BT, 'R_RBFSb','reactionName', 'Riboflavin synthase', 'reactionformula', '2 dmlz[C_c] ->	ribflv[C_c] + 4r5au[C_c] + h[C_c]','geneRule','( NP_813164_2 or ( NP_810230_1 and NP_813164_2 ))','lowerBound',0,'upperBound',1000);%original 2 dmlz[C_c] ->	ribflv[C_c] + 4r5au[C_c]

model_BT = removeRxns(model_BT, {'R_SALCHS4FER2'});
model_BT=addReaction(model_BT, 'R_SALCHS4FER2','reactionName', 'Salmochelin S4 Fe III reduction Fe III unloading', 'reactionformula', 'fmnh2[C_c] + 2 salchs4fe[C_c] -> 3 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 salchs4[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original fmnh2[C_c] + 2 salchs4fe[C_c] 	->	2 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 salchs4[C_c]

model_BT = removeRxns(model_BT, {'R_UDCPAT'});
model_BT=addReaction(model_BT, 'R_UDCPAT','reactionName', 'Abequosyl transferase', 'reactionformula', 'cdpabeq[C_c] + udcpgrm[C_c] -> 3 h[C_c] + cdp[C_c] + udcpo4[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original cdpabeq[C_c] + udcpgrm[C_c] 	->	h[C_c] + cdp[C_c] + udcpo4[C_c]

model_BT = removeRxns(model_BT, {'R_OOR3r'});
model_BT=addReaction(model_BT, 'R_OOR3r','reactionName', '2-oxoglutarate synthase (rev)', 'reactionformula', 'coa[C_c] + akg[C_c] + fdxox[C_c] -> h[C_c] + co2[C_c] + succoa[C_c] + fdxrd[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original h[C_c] + coa[C_c] + akg[C_c] + fdxox[C_c] 	->	co2[C_c] + succoa[C_c] + fdxrd[C_c]

model_BT = removeRxns(model_BT, {'R_ALCD2y'});
model_BT=addReaction(model_BT, 'R_ALCD2y','reactionName', 'ALCD2y', 'reactionformula', 'nad[C_c] + etoh[C_c]  <=> h[C_c] + nadh[C_c] + acald[C_c]','geneRule','','lowerBound',-1000,'upperBound',1000);%NP_813423_1

model_BT = removeRxns(model_BT,{'R_ACOAD1','R_ACOAD4_1','R_ACOAD2'});% These reactions used the wrong cofactor (NAD instead of FAD)

% Reactions modified because unlikely thermodynamically feasible (Nikita Report)
    
    model_BT = changeRxnBounds(model_BT,{'R_NDPK1';'R_NDPK2';'R_NDPK3';'R_NDPK4';'R_NDPK5';'R_NDPK6';'R_NDPK7';'R_NDPK8';'R_ADK4'},repelem(0,9),'l'); %Nucleoside-diphosphatekinase
    model_BT = changeRxnBounds(model_BT,'R_ADK1',0,'l');
    model_BT = removeRxns(model_BT, {'R_GK2','R_ADK2','R_ADK3','R_ALATA_L','R_DADK'});%,'R_ADK4',,'R_NDPK9'
    model_BT=addReaction(model_BT, 'R_GK2','reactionName', 'Guanylate kinase  GMPdATP', 'reactionformula', 'datp[C_c] + dgmp[C_c]  <=> dgdp[C_c] + dadp[C_c]','lowerBound',-1000,'upperBound',1000); %,'geneRule','NP_810922_1' https://enzyme.expasy.org/EC/2.7.4.8
    model_BT = changeRxnBounds(model_BT,{'R_GK2'},0,'b');
   
%% check for duplicated reactions (method 'FR')

 % Remove duplicated reactions and loop reactions

    model_BT = removeRxns(model_BT,{'R_ALAD_syn','R_CYSTAi','R_LDAPAT','R_FOMETRi','R_GLYAT','R_GLYCL_2','R_LGTHL','R_IG3PS_1','R_ICH','R_INDOLEt2pp','R_PGK_1','R_PGM_1','R_TDPDRR','R_XYLI2i','R_NTP1','R_IZPN','R_METS','R_PRAIS','R_NCTPPRT','R_BTS_1','R_NAt3_1'}); % NTP1 was atp[C_c] + h2o[C_c] => adp[C_c] + pi[C_c] + h[C_c];check with MetaCyc 
    model_BT = changeRxnBounds(model_BT, {'R_FERIRDe';'R_FRD2rpp';'R_GLYCL';'R_GLUDxi';'R_NAt3pp'},[-1000,-1000,-1000,-1000,-1000],'l');% 
    model_BT = changeRxnBounds(model_BT, {'R_FERIRDe';'R_FRD2rpp'},[0,1000],'u');% 
    model_BT=addReaction(model_BT, 'R_NADHDH2','reactionName', 'R_NADHDH2', 'reactionformula', '3 h[C_c] + nadh[C_c] + mqn8[C_c] -> nad[C_c] + 2 h[C_e] + mql8[C_c]','geneRule','','lowerBound',0,'upperBound',1000);
    model_BT = removeRxns(model_BT,{'R_NH4tpp_1','R_G6PI','R_G6PBDH','R_H2Ot','R_NH4t','R_GLCURtex','R_GLCURt2rpp','R_FUCtex','R_FUCtpp','R_GALURt2r','R_CO2t','R_ADEt2rpp','R_ADEtex','R_PIt2r','R_INSTtex','R_INOSTt4pp'}); %,'PROD2',,'R_HYDFDN','R_DGK1',,'R_INS2D','R_PROt4'
    model_BT = changeRxnBounds(model_BT,{'R_G3PD1';'R_G3PD2';'R_GLUDxi';'R_GLUDy';'R_HYDFDN2r'},[0,0,0,0,0],'u'); % https://biocyc.org/gene?orgid=BTHE&id=G13PU-7065;https://biocyc.org/gene?orgid=BTHE&id=G13PU-7065;https://biocyc.org/gene?orgid=BTHE&id=G13PU-6911-MONOMER#tab=RXNS;https://biocyc.org/gene?orgid=BTHE&id=G13PU-6911-MONOMER#tab=RXNS
    model_BT = changeRxnBounds(model_BT,{'R_ASPT';'R_DHDPRx_r';'R_THRAi';'R_DHDPRy';'R_ASPt2_3'},[-1000,0,-1000,0,0],'l'); %https://biocyc.org/gene?orgid=BTHE&id=G13PU-7716-MONOMER# BT_2755;https://biocyc.org/gene?orgid=ECOLI&id=DIHYDROPICRED-MONOMER (In Ecoli);https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-6563(In Ecoli);https://biocyc.org/gene?orgid=ECOLI&id=DIHYDROPICRED-MONOMER (In Ecoli)
    model_BT = changeRxnBounds(model_BT,{'R_ASPTA';'R_FRNDPR2r_1';'R_GLYCL'},[0,0,0],'b'); % https://biocyc.org/gene?orgid=BTHE&id=G13PU-9292-MONOMER#tab=RXNS

%% Add reactions for ethanol production

model_BT=addReaction(model_BT, 'R_ADLH','reactionName', 'acetaldehyde:NAD+ oxidoreductase (CoA-acetylating)', 'reactionformula', 'acald[C_c] + coa[C_c] + nad[C_c] -> nadh[C_c] + accoa[C_c] + h[C_c]','lowerBound',-1000,'upperBound',1000);%original 'geneRule','BT_3115'

%% Change the ATPM name for 
model_BT = removeRxns(model_BT, {'R_ATPM'});
model_BT=addReaction(model_BT, 'R_ngam','reactionName', 'Non-growth ATP maintenance', 'reactionformula', 'atp[C_c] + h2o[C_c]  -> adp[C_c] + pi[C_c] + h[C_c]','geneRule','spontaneous','lowerBound',0,'upperBound',1000);

%% Added to the model
  % Add metabolites 
  met={'3mb[C_c]','acgal[C_e]','acgal[C_c]','acgal6P[C_c]','gal6P[C_c]','acgam[C_e]','cys__L[C_e]','gly[C_e]','his__L[C_e]','isobuta[C_c]','isobuta[C_e]','lys__L[C_e]','met__L[C_e]','phe__L[C_e]','ser__L[C_e]','thr__L[C_e]','trp__L[C_e]','tyr__L[C_e]','btn[C_e]','s[C_e]','pydam[C_e]'};%,'hco3[C_e]','pnto__R[C_e]','acgal[C_c]',
        
        model_BT= addMultipleMetabolites(model_BT,met);
    
 %% Add formulas and Names of metabolites added before
    
             % Add Exchange reactions
        rxnToAdd.metList     = {'acgal[C_e]','acgam[C_e]','cys__L[C_e]','gly[C_e]','his__L[C_e]','isobuta[C_e]','lys__L[C_e]','met__L[C_e]','phe__L[C_e]','ser__L[C_e]','thr__L[C_e]','trp__L[C_e]','tyr__L[C_e]','btn[C_e]','s[C_e]','pydam[C_e]'};%'na1[C_e]','pnto__R[C_e]','hco3[C_e]'
        rxnToAdd.lb    = repelem(-1000,length(rxnToAdd.metList));
        rxnToAdd.ub    = repelem(1000,length(rxnToAdd.metList));
        model_BT = addExchangeRxn_v2(model_BT,rxnToAdd.metList,rxnToAdd.lb,rxnToAdd.ub); %the function has been modified to add 'R_' in front of exchange reaction names
    
            % Add transport reactions 
                
                % Based on Thiele I, et al (2017)
                model_BT=addReaction(model_BT, 'R_Cystr','reactionName', 'L-Cysteine reversible transport via proton symport', 'reactionformula', 'h[C_e] + cys__L[C_e] <=> h[C_c] + cys__L[C_c]'); 
                model_BT=addReaction(model_BT, 'R_Glyt2r','reactionName', 'Glycine reversible transport via proton symport', 'reactionformula', '2 h[C_e] + gly[C_e] <=> 2 h[C_c] + gly[C_c]'); 
                model_BT=addReaction(model_BT, 'R_Argt2r','reactionName', 'Arginine reversible transport via proton symport', 'reactionformula', 'h[C_e] + arg__L[C_e] <=> h[C_c] + arg__L[C_c]'); 
                model_BT=addReaction(model_BT, 'R_Hist2r','reactionName', 'L-Histidine reversible transport via proton symport', 'reactionformula', '2 h[C_e] + his__L[C_e] <=> 2 h[C_c] + his__L[C_c]');
                model_BT=addReaction(model_BT, 'R_Lyst2r','reactionName', 'L-Lysine reversible transport via proton symport', 'reactionformula', '2 h[C_e] + lys__L[C_e] <=> 2 h[C_c] + lys__L[C_c]');
                model_BT=addReaction(model_BT, 'R_Met2r','reactionName', 'L-Methionine reversible transport via proton symport', 'reactionformula', 'h[C_e] + met__L[C_e] <=> h[C_c] + met__L[C_c]');
                model_BT=addReaction(model_BT, 'R_PHEt2r','reactionName', 'L-Phenylalanine reversible transport via proton symport', 'reactionformula', 'h[C_e] + phe__L[C_e] <=> h[C_c] + phe__L[C_c]'); 
                model_BT=addReaction(model_BT, 'R_SERt2r','reactionName', 'L-Serine reversible transport via proton symport', 'reactionformula', 'h[C_e] + ser__L[C_e] <=> h[C_c] + ser__L[C_c]'); 
                model_BT=addReaction(model_BT, 'R_THRt2r','reactionName', 'L-Threonine reversible transport via proton symport', 'reactionformula', 'h[C_e] + thr__L[C_e] <=> h[C_c] + thr__L[C_c]'); 
                model_BT=addReaction(model_BT, 'R_TRPt2r','reactionName', 'L-Tryptophan reversible transport via proton symport', 'reactionformula', 'h[C_e] + trp__L[C_e] <=> h[C_c] + trp__L[C_c]');
                model_BT=addReaction(model_BT, 'R_Leut2r','reactionName', 'L-Leucine reversible transport via proton symport', 'reactionformula', 'h[C_e] + leu__L[C_e] <=> h[C_c] + leu__L[C_c]'); 
                model_BT=addReaction(model_BT, 'R_Valt2r','reactionName', 'L-Valine reversible transport via proton symport', 'reactionformula', 'h[C_e] + val__L[C_e] <=> h[C_c] + val__L[C_c]'); 
                model_BT=addReaction(model_BT, 'R_Tyrt2r','reactionName', 'L-Tyrosine reversible transport via proton symport', 'reactionformula', 'h[C_e] + tyr__L[C_e] <=> h[C_c] + tyr__L[C_c]');            
                model_BT=addReaction(model_BT, 'R_Isobut','reactionName', 'Isobutyrate transport', 'reactionformula', 'isobuta[C_e] + h[C_e] <=> isobuta[C_c] + h[C_c]'); % Isobutyric acid transporter => based on butyrate transporter
                model_BT = removeRxns(model_BT, {'R_PPAt4pp'});
                model_BT=addReaction(model_BT, 'R_PPAt4pp','reactionName', 'Propionate transport', 'reactionformula', 'na1[C_p] + ppa[C_p] <=> na1[C_c] + ppa[C_c]','geneRule','','lowerBound',-1000,'upperBound',1000); % originally monodirectional in the model
                model_BT=addReaction(model_BT, 'R_cittr','reactionName', 'Citrate transport', 'reactionformula', 'cit[C_e] + 3 na1[C_e] <=> cit[C_c] + 3 na1[C_c]'); % 100% ID; 100% Query cover/ tr|A0A0P0FBG3|A0A0P0FBG3_BACT4 Citrate transporter OS=Bacteroides thetaiotaomicron OX=818 GN=BSIG_0234 PE=4 SV=1; Transport mechanism not described on BT. Assume to be symport                 
                model_BT=addReaction(model_BT, 'R_ILEt2r','reactionName', 'Isoleucine reversible transport via proton symport', 'reactionformula', 'h[C_p] + ile__L[C_p] <=> h[C_c] + ile__L[C_c]');             
                model_BT=addReaction(model_BT, 'R_btnt','reactionName', 'Biotin transport', 'reactionformula', 'atp[C_c] + btn[C_e] + h2o[C_c]  -> adp[C_c] + btn[C_c] + h[C_c] + pi[C_c]'); 
                model_BT=addReaction(model_BT, 'R_St','reactionName', 'sulfur transport', 'reactionformula', 's[C_e] -> s[C_c]'); 
                model_BT=addReaction(model_BT, 'R_Ibutsyn','reactionName', 'Isobutyrate synthesis', 'reactionformula', '4mop[C_c] + h2o[C_c] -> isobuta[C_c] + ac[C_c]'); %Free gibbs = -14.957 kcal/mol. Described in Firmicutes/Clostridium Sporogens. (METACYC leucine degradation II).   
                model_BT=addReaction(model_BT, 'R_PDXt','reactionName', 'Pyridoxamine transport', 'reactionformula', 'pydam[C_e] <=> pydam[C_c]');  
                model_BT=addReaction(model_BT, 'R_galtr','reactionName', 'Galactose transport', 'reactionformula', 'atp[C_c] + gal[C_p] + h2o[C_c]  -> adp[C_c] + gal[C_c] + h[C_c] + pi[C_c]','lowerBound',0,'upperBound',1000);% ,'geneRule','NP_813222_1?' %uniprot.org/uniprot/Q89ZR5
                model_BT=addReaction(model_BT, 'R_PyrDH','reactionName', 'Pyruvate dehydrogenase', 'reactionformula', 'coa[C_c] + pyr[C_c] + nad[C_c]  -> accoa[C_c] + co2[C_c] + nadh[C_c]','lowerBound',0,'upperBound',1000);% ,'geneRule','BT_1820',https://biocyc.org/gene?orgid=BTHE&id=G13PU-6757-MONOMER

           %Add N-Acetylglucosamine N-acetylglucosamine permease MFS-type transporter but uncharacterized https://pubseed.theseed.org/?page=Annotation&feature=fig|226186.12.peg.4733 and https://prosite.expasy.org/PDOC50850
           
            model_BT=addReaction(model_BT, 'R_NAGp','reactionName', 'Permease', 'reactionformula', 'h[C_e] + acgam[C_e] -> h[C_c] + acgam[C_c]','lowerBound',0,'upperBound',1000); 
            model_BT=addReaction(model_BT, 'R_PTSI_AgaVWEF','reactionName', 'NAcGal transport', 'reactionformula', 'acgal[C_e] + h[C_e] <=> acgal[C_c] + h[C_c]','lowerBound',0,'upperBound',1000); %
            model_BT=addReaction(model_BT, 'R_AcgalP','reactionName', 'Acgal kinase', 'reactionformula', 'acgal[C_c] + atp[C_c]  -> acgal6P[C_c] + adp[C_c] + h[C_c]','lowerBound',0,'upperBound',1000); %
            model_BT=addReaction(model_BT, 'R_AgaAI','reactionName', 'Deacetylase', 'reactionformula', 'acgal6P[C_c] + h2o[C_c] -> gal6P[C_c] + ac[C_c]','lowerBound',0,'upperBound',1000); %
            model_BT=addReaction(model_BT, 'R_AgaS','reactionName', 'D-galactosamine-6-phosphate deaminase/isomerase', 'reactionformula', 'gal6P[C_c] + h2o[C_c] -> tag6p__D[C_c] + nh4[C_c]','lowerBound',0,'upperBound',1000); %
            model_BT=addReaction(model_BT, 'R_Succ_Asp_symport','reactionName', 'Succinate transport', 'reactionformula', 'succ[C_c] + asp__L[C_e] <=> asp__L[C_c] + succ[C_e]','lowerBound',-1000,'upperBound',1000); %  W/O fumarate prediction of GR was unfeasible fumarate symporter => change into aspartate transporter; geneRule: BT_2756
            model_BT=addReaction(model_BT, 'R_Isoval','reactionName', 'IsovalerateCoA hydrolysis', 'reactionformula', 'ivcoa[C_c] + h2o[C_c] -> coa[C_c] + h[C_c] + 3mb[C_c]'); %https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-8185
            model_BT=addReaction(model_BT, 'R_Isovaltr','reactionName', 'IsovalerateCoA diffusion periplasmic', 'reactionformula', '3mb[C_c] <=> 3mb[C_p]'); %
    
             %Two cytochrome C metabolites were in the model   
             model_BT=removeMetabolites(model_BT, {'focytC[C_c]'})
             model_BT = addMetabolite(model_BT,'focytred[C_c]','metFormula','C42FeH54N8O6S2','Charge',0,'metName','Ferrocytochrome C')

             %'metNotes','BioCyc: META:Cytochromes-C-Reduced;SEED Compound: cpd00110;UniPathway Compound: UPC00126;KEGG Compound: C00126;MetaNetX (MNX) Chemical: MNXM537');

             model_BT = removeRxns(model_BT, {'R_CYOO2pp','R_CYTBCYTC'});
             model_BT=addReaction(model_BT, 'R_CYOO2pp','reactionName', 'Cytochrome-c oxidase (2 protons translocated) periplasm', 'reactionformula', '0.5 o2[C_c] + 4 h[C_c] + 2 focytred[C_c]	->	h2o[C_c] + 2 h[C_p] + 2 ficytC[C_c]','lowerBound',0,'upperBound',1000); %
             model_BT=addReaction(model_BT, 'R_CYTBCYTC','reactionName', 'Cytochrome b Cytochrome c electron transfer', 'reactionformula', 'focytc[C_c] + ficytC[C_c] 	<=>	focytred[C_c] + ficytb[C_c]','lowerBound',-1000,'upperBound',1000); %

%% Add formulas and Names of metabolites added before

        Formulas= {'C5H9O2','C8H15NO6','C8H15NO6','C8H14NO9P','C6H13NO8P','C8H15NO6','C3H7NO2S','C2H5NO2','C6H9N3O2','C4H8O2','C4H8O2','C6H15N2O2','C5H11NO2S','C9H11NO2','C3H7NO3','C4H9NO3','C11H12N2O2','C9H11NO3','C10H15N2O3S','S','C8H13N2O2'};%,'CHO3','C9H15NO5','C8H15NO6',
        Names= {'3-Methylbutanoic acid','N-Acetyl-D-galactosamine','N-Acetyl-D-galactosamine','N-Acetyl-D-galactosamine6P','D-galactosamine-6-phosphate','N-Acetyl-D-glucosamine','L-Cysteine','Glycine','L-Histidine','Isobutyrate','Isobutyrate','L-Lysine','L-Methionine','L-Phenylalanine','L-Serine','L-Threonine','L-Tryptophan','L-Tyrosine','Biotin','sulfur','Pyridoxamine'};%,'Bicarbonate','Pantothenate','N-Acetyl-D-galactosamine',

        for i=1:numel(met)            
            pos = strcmp(model_BT.mets,char(met(i)));
            k= find(pos);
            model_BT.metNames(k)= Names(i);
            model_BT.metFormulas(k) = Formulas(i);
        end

%% Model constrains

            %Amino Acids
            AA= {'R_EX_ala__L_e';'R_EX_arg__L_e';'R_EX_asn__L_e';'R_EX_asp__L_e';'R_EX_citr__L_e';'R_EX_cys__L_e';'R_EX_gln__L_e';'R_EX_glu__L_e';'R_EX_gly_e';'R_EX_his__L_e';'R_EX_ile__L_e';'R_EX_leu__L_e';'R_EX_lys__L_e';'R_EX_met__L_e';'R_EX_phe__L_e';'R_EX_pro__L_e';'R_EX_ser__L_e';'R_EX_thr__L_e';'R_EX_trp__L_e';'R_EX_tyr__L_e';'R_EX_val__L_e'};
            model_BT = changeRxnBounds(model_BT,AA,[repelem(0,length(AA))],'u'); %Hypothesis: BT is not producing any Amino acids
            model_BT = changeRxnBounds(model_BT,AA,[repelem(-1000,length(AA))],'l'); 
            
            AA_other={'R_EX_ala__D_e';'R_EX_orn_e';'R_EX_progly_e'};
            model_BT = changeRxnBounds(model_BT,AA_other,[repelem(0,length(AA_other))],'u'); 
            model_BT = changeRxnBounds(model_BT,AA_other,[repelem(0,length(AA_other))],'l'); 
           
            % Nucleic Acids
            NA= {'R_EX_adn_e';'R_EX_gsn_e';'R_EX_uri_e';'R_EX_thymd_e'};
            model_BT = changeRxnBounds(model_BT,NA,[repelem(0,length(NA))],'u'); 
            model_BT = changeRxnBounds(model_BT,NA,[repelem(0,length(NA))],'l');
                        
            NA_deo= {'R_EX_dad_2_e';'R_EX_dcyt_e';'R_EX_dgsn_e';'R_EX_duri_e'};
            model_BT = changeRxnBounds(model_BT,NA_deo,[repelem(0,length(NA_deo))],'u'); 
            model_BT = changeRxnBounds(model_BT,NA_deo,[repelem(-1000,length(NA_deo))],'l'); % Hypothesis: RNA are more abondant and less stable in cell than DNA making nucleoside more available.
            
            NA_other= {'R_EX_3amp_e';'R_EX_3cmp_e';'R_EX_3gmp_e';'R_EX_ade_e';'R_EX_cytd_e';'R_EX_din_e';'R_EX_gua_e';'R_EX_ins_e';'R_EX_thym_e';'R_EX_ura_e';'R_EX_xtsn_e'};
            model_BT = changeRxnBounds(model_BT,NA_other,[repelem(0,length(NA_other))],'u'); 
            model_BT = changeRxnBounds(model_BT,NA_other,[repelem(0,length(NA_other))],'l');
            
            % Carbohydrates
            CH= {'R_EX_abt__L_e';'R_EX_amylose300_e';'R_EX_Larab_e';'R_EX_araban__L_e';'R_EX_arab__D_e';'R_EX_arab__L_e';'R_EX_cell6_e';'R_EX_dha_e';'R_EX_drib_e';'R_EX_galt_e';'R_EX_galur_e';'R_EX_glc__D_e';'R_EX_glyald_e';'R_EX_lcts_e';'R_EX_lyx__L_e';'R_EX_pullulan1200_e';'R_EX_raffin_e';'R_EX_rib__D_e';'R_EX_rmn_e';'R_EX_xyl__D_e';'R_EX_xyl3_e';'R_EX_xylan4_e';'R_EX_xylb_e'};
            model_BT = changeRxnBounds(model_BT,CH,[repelem(0,length(CH))],'b');
          
            % Mucin Component (carbohydrate part)
            MCP= {'R_EX_fuc__L_e';'R_EX_gal_e';'R_EX_acgam_e';'R_EX_acgal_e';'R_EX_anhgm_e';'R_EX_gal_bD_e'}; %According to Thiele mannose is not in mucin 
            model_BT = changeRxnBounds(model_BT,MCP,[repelem(-1000,4),0,0],'l');
            model_BT = changeRxnBounds(model_BT,MCP,[repelem(0,length(MCP))],'u');

            %Ions and Vitamins
            IV={'R_EX_btn_e';'R_EX_cobalt2_e';'R_EX_pydam_e';'R_EX_s_e';'R_EX_ca2_e';'R_EX_cl_e';'R_EX_cu2_e';'R_EX_fe2_e';'R_EX_fe3_e';'R_EX_h_e';'R_EX_k_e';'R_EX_mg2_e';'R_EX_mn2_e';'R_EX_nh4_e';'R_EX_no2_e';'R_EX_no3_e';'R_EX_pi_e';'R_EX_so4_e';'R_EX_zn2_e'};
            model_BT = changeRxnBounds(model_BT,IV,[repelem(1000,length(IV))],'u'); 
            model_BT = changeRxnBounds(model_BT,IV,[repelem(-1000,length(IV))],'l'); 
            
            %model_BT = changeRxnBounds(model_BT,{'R_EX_pnto__R_e'},0,'u'); 
            model_BT = changeRxnBounds(model_BT, 'R_EX_no3_e',0,'l');

            % Gaz and water
            Gaz={'R_EX_co2_e';'R_EX_h2s_e';'R_EX_o2_e';'R_EX_h2_e';'R_EX_h2o_e'};
            model_BT = changeRxnBounds(model_BT,Gaz,[1000,0,repelem(0,2),1000],'u'); %(Predicted growth rate 0.1428 /h if co2 UB=7)
            model_BT = changeRxnBounds(model_BT,Gaz,[0,0,repelem(0,2),-1000],'l');

            % Lipid and Membrane components
            LipMc={'R_EX_chol_e';'R_EX_chols_e';'R_EX_dca_e';'R_EX_enlipa_e';'R_EX_g3pe_e';'R_EX_g3pg_e';'R_EX_glcur_e';'R_EX_glucan4_e';'R_EX_glucan6_e';'R_EX_hxan_e';'R_EX_inost_e';'R_EX_LalaDglu_e';'R_EX_kdo2lipid4_e';'R_EX_LalaLglu_e';'R_EX_LalaDgluMdap_e';'R_EX_LalaDgluMdapDala_e';'R_EX_malt_e';'R_EX_malttr_e';'R_EX_manttr_e';'R_EX_octa_e';'R_EX_udcpdp_e';'R_EX_udcpo5_e';'R_EX_udcpp_e'};
            model_BT = changeRxnBounds(model_BT,LipMc,[repelem(0,length(LipMc))],'b');
                                   
            % Siderophore (enterobactin); antioxidant (glutathion)
            Sid={'R_EX_enter_e';'R_EX_feenter_e';'R_EX_salchs4_e';'R_EX_salchs4fe_e';'R_EX_gthox_e';'R_EX_gthrd_e'};
            model_BT = changeRxnBounds(model_BT,Sid,[repelem(0,length(Sid)-2),-1000,-1000],'l');
            model_BT = changeRxnBounds(model_BT,Sid,[repelem(0,length(Sid)-2),1000,1000],'u');

            % Unclassified
            model_BT = changeRxnBounds(model_BT,{'R_EX_acald_e'},0,'u'); %acald and fumarate are in the media
            model_BT = changeRxnBounds(model_BT,{'R_EX_acald_e'},0,'l'); 
           
            model_BT = changeRxnBounds(model_BT,{'R_EX_fum_e'},0,'u'); %acald and fumarate are in the media
            model_BT = changeRxnBounds(model_BT,{'R_EX_fum_e'},0,'l'); 
           
            Unclass={'R_EX_id3acald_e';'R_EX_12ppd__S_e';'R_EX_14glucan_e';'R_EX_2pglyc_e';'R_EX_4abut_e';'R_EX_4hbz_e';'R_EX_4hphac_e';'R_EX_alltn_e';'R_EX_aso3_e';'R_EX_aso4_e';'R_EX_cgly_e';'R_EX_etha_e';'R_EX_glyb_e';'R_EX_glyc_e';'R_EX_glyc3p_e';'R_EX_glycogen1500_e';'R_EX_h2o2_e';'R_EX_lmn2_e';'R_EX_lmn30_e';'R_EX_meoh_e';'R_EX_mso3_e';'R_EX_n2o_e';'R_EX_no_e';'R_EX_pac_e';'R_EX_pyr_e';'R_EX_spmd_e';'R_EX_starch1200_e';'R_EX_tartr__D_e';'R_EX_tol_e';'R_EX_xan_e'};
            model_BT = changeRxnBounds(model_BT,Unclass,[repelem(0,length(Unclass))],'b');
          
            % Required fluxes to optimize Biomass
            model_BT = changeRxnBounds(model_BT,'R_EX_glyclt_e',1,'u'); %
            model_BT = changeRxnBounds(model_BT,'R_EX_glyclt_e',0,'l'); 

            % Experimental rates

            model_BT = changeRxnBounds(model_BT, 'R_EX_ac_e',4.13,'u');%
            model_BT = changeRxnBounds(model_BT, 'R_EX_ac_e',3.34,'l');

            model_BT = changeRxnBounds(model_BT, 'R_EX_but_e',0,'u');
            model_BT = changeRxnBounds(model_BT, 'R_EX_but_e',0,'l');

            model_BT = changeRxnBounds(model_BT,'R_EX_cit_e',0.37,'u');% Original model can only import citrate (ATP dpt transpoter). Transport reaction has been added
            model_BT = changeRxnBounds(model_BT,'R_EX_cit_e',0,'l');

            model_BT = changeRxnBounds(model_BT, 'R_EX_succ_e',1.54,'u');
            model_BT = changeRxnBounds(model_BT, 'R_EX_succ_e',1.30,'l');

            model_BT = changeRxnBounds(model_BT, 'R_EX_for_e',0.54,'u');
            model_BT = changeRxnBounds(model_BT, 'R_EX_for_e',0.39,'l');

            model_BT = changeRxnBounds(model_BT, 'R_EX_etoh_e',3.99,'u');
            model_BT = changeRxnBounds(model_BT, 'R_EX_etoh_e',3.3,'l');
            
            model_BT = changeRxnBounds(model_BT, 'R_EX_ppa_e',3.87,'u');%   
            model_BT = changeRxnBounds(model_BT, 'R_EX_ppa_e',3.16,'l');%

            model_BT = changeRxnBounds(model_BT, 'R_EX_lac__D_e',-0.80,'u');% 
            model_BT = changeRxnBounds(model_BT, 'R_EX_lac__D_e',-0.92,'l');% 

            model_BT = changeRxnBounds(model_BT, 'R_EX_isobuta_e',0.21,'u');% 
            model_BT = changeRxnBounds(model_BT, 'R_EX_isobuta_e',0.16,'l');%         
            
            model_BT = changeRxnBounds(model_BT, 'R_EX_mal__L_e',-0.17,'u');% 
            model_BT = changeRxnBounds(model_BT, 'R_EX_mal__L_e',-0.22,'l');% 
            
            model_BT = changeRxnBounds(model_BT, 'R_EX_3mb_e',0.79,'u');% isovalerate
            model_BT = changeRxnBounds(model_BT, 'R_EX_3mb_e',0.70,'l');% 
                
%% Biomass reaction modification.
     % Biomass GAM modified to reach experimental growth rate
     model_BT = addReaction(model_BT,'Growth','reactionName', 'Biomass reaction','reactionFormula','0.000215957 10fthf[C_c] + 26 atp[C_c] + 0.1 amp[C_c] + 0.254849 glu__L[C_c] + 26 h2o[C_c] + 0.00177318 nad[C_c] + 0.000557809 coa[C_c] + 0.000432883 nadp[C_c] + 0.124 gmp[C_c] + 0.000215957 amet[C_c] + 0.436303 leu__L[C_c] + 0.233442 asp__L[C_c] + 0.000215957 fad[C_c] + 0.497466 ala__L[C_c] + 0.00504062 cl[C_c] + 0.000215957 ribflv[C_c] + 0.254849 gln__L[C_c] + 0.0886878 cys__L[C_c] + 0.208977 ser__L[C_c] + 0.245675 thr__L[C_c] + 0.0917461 his__L[C_c] + 0.000215957 thmpp[C_c] + 0.074 cmp[C_c] + 0.0136 damp[C_c] + 0.0143 dgmp[C_c] + 0.59329 gly[C_c] + 0.148832 met__L[C_c] + 0.000215957 thf[C_c] + 0.233442 asn__L[C_c] + 0.000215957 pydx5p[C_c] + 0.286451 arg__L[C_c] + 9.68419e-05 mql8[C_c] + 0.214074 pro__L[C_c] + 0.332324 lys__L[C_c] + 0.00650293 fe2[C_c] + 0.179414 phe__L[C_c] + 0.004201 so4[C_c] + 0.189029 k[C_c] + 0.281354 ile__L[C_c] + 0.000686609 cu2[C_c] + 0.00504062 ca2[C_c] + 0.00840103 mg2[C_c] + 0.000669178 mn2[C_c] + 9.68419e-05 cobalt2[C_c] + 0.000330231 zn2[C_c] + 0.00756142 fe3[C_c] + 0.0802 ump[C_c] + 0.0137 dcmp[C_c] + 0.000215957 mlthf[C_c] + 0.0135 dtmp[C_c] + 0.0550478 trp__L[C_c] + 0.409798 val__L[C_c] + 0.133541 tyr__L[C_c] + 0.0968419 uaagmda[C_c] 	->	26 adp[C_c] + 26 pi[C_c] + 26 h[C_c]');

    %% Optimization for biomass

    model_BT=changeObjective(model_BT,'Growth'); %Predicted is 0.1224
    FBAsolutionBT=optimizeCbModel(model_BT,'max','one',false);
    printFluxVector(model_BT, FBAsolutionBT.x, 'true');

%% Export model as SMBL containing constrains

    writeCbModel(model_BT,'sbml','[your_path]/Bacteroides_thetaiotaomicron_VPI_5482',{'C_c','C_p','C_e'},{'cytosol','periplasm','extracellular space'}) %SBML level = 2 by default, sbmlVersion = 1 by default

    %% Random Sampling (RAVEN)
        %not working so well with COBRATOOLBOX

        cd '[your_path]';
        model_BTrav=importModel('Bacteroides_thetaiotaomicron_VPI_5482.xml')%1error

    % Fix measured exchange fluxes around 10% of the value from FBA
    exIdx = getIndexes(model_BTrav,{'EX_ac_e','EX_cit_e','EX_succ_e','EX_ppa_e','EX_lac__D_e','EX_for_e','EX_etoh_e','EX_mal__L_e', 'EX_but_e','EX_isobuta_e','EX_3mb_e','Growth','ngam'},'rxns');%

    fluxes = sol.x(exIdx);
    model_BTrav = setParam(model_BTrav,'var',exIdx,fluxes,10);    
    [~, goodrxn] = randomSampling(model_BTrav,1,true,true,true); 
    rs = randomSampling(model_BTrav,5000,true,true,true,goodrxn,false); % 5000 Sampling generated
    fluxMean = full(mean(rs,2));
    fluxSD = full(std(rs,0,2));
    
    %% Save FBA predicted fluxes
    
    cd '[your_path]';
 
    out=table(model_BTrav.rxns,model_BTrav.rxnNames,constructEquations(model_BTrav),sol.x);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Flux'};
    out(out.Flux==0,:) =[];
    writetable(out,'BT_FBA.xlsx');


    %%  Save FVA predicted fluxes (mean+SD) and random sampling matrix
    clear out;out = fluxMean;
    clear out2;out2 = fluxSD;
    
    cd '[your_path]';

    out =table(model_BTrav.rxns,model_BTrav.rxnNames, constructEquations(model_BTrav), out ,out2);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Mean_RS' 'SD_RS'};
    out(out.Mean_RS==0,:) =[];
    writetable(out,'BT_RS.xlsx');

    rs(~any(rs,2),:) =[];
    RS_df=   array2table(full(rs));
    RS_df = [out(:,[1:3]),RS_df];
    writetable(RS_df,'BT_RS_5000samples.xlsx');

%% NADH and NADPH reactions
i= 'Nicotinamide adenine dinucleotide phosphate - reduced';
    [fluxes, rxnIdx] = getMetProduction(model_BTrav,i,fluxMean,true);  
    [SD, rxnIdx] = getMetProduction(model_BTrav,i,fluxSD,true);   
    clear out3;
    out3.rxns    = model_BTrav.rxns(rxnIdx);
    out3.rxnNames= model_BTrav.rxnNames(rxnIdx);
    out3.rxnEqns = constructEquations(model_BTrav,rxnIdx);
    out3.fluxes  = num2cell(fluxes);
    out3.SD  = num2cell(SD); 
    out3 = table(out3.rxns,out3.rxnNames,out3.rxnEqns,out3.fluxes,out3.SD);
    out3.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'RS_Flux' 'RS_SD'};
    writetable(out3,'BT_NADPH_RS.xlsx');
    
i= 'Nicotinamide adenine dinucleotide - reduced';
    [fluxes, rxnIdx] = getMetProduction(model_BTrav,i,fluxMean,true);
    [SD, rxnIdx] = getMetProduction(model_BTrav,i,fluxSD,true);   
    clear out4;
    out4.rxns    =model_BTrav.rxns(rxnIdx);
    out4.rxnNames= model_BTrav.rxnNames(rxnIdx);
    out4.rxnEqns = constructEquations(model_BTrav,rxnIdx);
    out4.fluxes  = num2cell(fluxes); 
    out4.SD  = num2cell(SD); 
    out4 = table(out4.rxns,out4.rxnNames,out4.rxnEqns,out4.fluxes,out4.SD);
    out4.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'RS_Flux' 'RS_SD'};
    writetable(out4,'BT_NADH_RS.xlsx');

    %% Individual depletion analysis
    EX_RS= out(contains(out.(1),'EX_'),:);%55 rows
    temp= EX_RS(EX_RS.Mean_RS <0,:); %44
    
    IV={'R_EX_btn_e';'R_EX_cobalt2_e';'R_EX_pydam_e';'R_EX_s_e';'R_EX_ca2_e';'R_EX_cl_e';'R_EX_cu2_e';'R_EX_fe2_e';'R_EX_fe3_e';'R_EX_h_e';'R_EX_k_e';'R_EX_mg2_e';'R_EX_mn2_e';'R_EX_nh4_e';'R_EX_no2_e';'R_EX_no3_e';'R_EX_pi_e';'R_EX_so4_e';'R_EX_zn2_e'};
     
    %Remove ions
    IV2=strrep(IV,'R_','');
    idx = ismember(temp.Reaction,IV2); %get indexes of ions and vitamins metabolites
    temp(idx,:)=[];%30
    
    % Remove exp metabolites
    Exp_meta={'EX_ac_e','EX_but_e','EX_cit_e','EX_succ_e','EX_for_e','EX_etoh_e','EX_ppa_e','EX_lac__D_e','EX_isobuta_e','EX_mal__L_e','EX_3mb_e'};
    idx2 = ismember(temp.Reaction,Exp_meta); %get indexes of ions and vitamins metabolites
    temp(idx2,:)=[];%28
    
    % Remove water
    temp(strcmp(temp.Reaction,'EX_h2o_e'),:)=[];%27
    
    temp.Reaction = strcat('R_', temp.Reaction);
     
    temp.name = strrep(temp.Reaction,'R_EX_','');

    cd '[your_path]';
    %first carbon source was threonine
    GR_prd_BT = table();
    
    for i=1:length(temp.Reaction)
        model_BT = changeRxnBounds(model_BT, temp.Reaction(i),0,'l');% 
        FBAsolution_temp=optimizeCbModel(model_BT,'max','one',false); % Optimize for GR        
        disp([temp.Reaction(i) FBAsolution_temp.f])
        GR_prd_BT.GR(i)= FBAsolution_temp.f;
        GR_prd_BT.Depelted_Metabolite(i)= temp.Reaction(i);
            if isnan(FBAsolution_temp.f)~=1
                outCOB=table(model_BT.rxns,model_BT.rxnNames,FBAsolution_temp.x);
                outCOB.Properties.VariableNames = {'Reaction' 'Description' 'Flux'};
                outCOB(outCOB.Flux==0,:) =[];
                 writetable(outCOB,string(strcat('BT_FBA_WO_',temp.Reaction(i),'.xlsx')));
            end
        model_BT = changeRxnBounds(model_BT,  temp.Reaction(i),-1000,'l');%
     end
 writetable(GR_prd_BT,string(strcat('BT_FBA_GR_indDepletion.xlsx')));

