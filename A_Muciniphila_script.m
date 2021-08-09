%initialization CobraToolbox and GUROBI as solver
initCobraToolbox
changeCobraSolver('gurobi','all')

%% Load the model

getFullPath(fullfile(cd, 'Akkermansia_muciniphila_ATCC_BAA_835.xml'));
cd '[your_path]';
model=readCbModel('Akkermansia_muciniphila_ATCC_BAA_835','fileType','SBML');

%% Rename the model
model_AM = model;

%% Curation
%Update metabolites formula

pos1 = strcmp(model_AM.metNames,'L-Cystathionine');
k= find(pos1);
model_AM.metFormulas(k);
model_AM.metFormulas(k)= {'C7H14N2O4S'};

pos2 = strcmp(model_AM.metNames,'L-Glutamyl-tRNA(Glu)');
k2= find(pos2);
model_AM.metFormulas(k2);
model_AM.metFormulas(k2)= {'C20H28N6O13PR'};

pos4 = strcmp(model_AM.metNames,'Pyridoxal 5''-phosphate');
k4= find(pos4);
model_AM.metFormulas(k4);
model_AM.metFormulas(k4)= {'C8H8NO6P'};

pos5 = strcmp(model_AM.metNames,'Folate');
k5= find(pos5);
model_AM.metFormulas(k5);
model_AM.metFormulas(k5)= {'C19H18N7O6'};

pos6 = strcmp(model_AM.metNames,'Sulfite');
k6= find(pos6);
model_AM.metFormulas(k6);
model_AM.metFormulas(k6)= {'HSO3'};

pos7 = strcmp(model_AM.metNames,'Oxidized thioredoxin');
k7= find(pos7);
model_AM.metFormulas(k7);
model_AM.metFormulas(k7)= {'C6H7NO2S2R2'};

pos8 = strcmp(model_AM.metNames,'Reduced thioredoxin');
k8= find(pos8);
model_AM.metFormulas(k8);
model_AM.metFormulas(k8)= {'C6H9NO2S2R2'};

pos9 = strcmp(model_AM.metNames,'Riboflavin');
k9= find(pos9);
model_AM.metFormulas(k9);
model_AM.metFormulas(k9)= {'C17H20N4O6'};

pos10 = strcmp(model_AM.metNames,'5-Methyltetrahydrofolate');
k10= find(pos10);
model_AM.metFormulas(k10);
model_AM.metFormulas(k10)= {'C20H24N7O6'};

% Update existing reactions
model_AM = removeRxns(model_AM, 'R_GLUDy');
model_AM=addReaction(model_AM, 'R_GLUDy','reactionName', 'Glutamate dehydrogenase (NADP)', 'reactionformula', 'nh4[C_c] + h[C_c] + nadph[C_c] + akg[C_c]  -> glu__L[C_c] + h2o[C_c] + nadp[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012421075_1');%original glu__L[C_c] + h2o[C_c] + nadp[C_c]  <=> nh4[C_c] + h[C_c] + nadph[C_c] + akg[C_c]

model_AM = removeRxns(model_AM, 'R_RBFK');
model_AM=addReaction(model_AM, 'R_RBFK','reactionName', 'Riboflavin kinase', 'reactionformula', 'atp[C_c] + ribflv[C_c]  -> adp[C_c] + 2 h[C_c] + fmn[C_c] ','lowerBound',0,'upperBound',1000,'geneRule','WP_012419799_1');%original atp[C_c] + ribflv[C_c]  -> adp[C_c] + h[C_c] + fmn[C_c] ; WP_012419799_1

model_AM = removeRxns(model_AM, 'R_APSR');
model_AM=addReaction(model_AM, 'R_APSR','reactionName', 'APSR', 'reactionformula', 'aps[C_c] + trdrd[C_c]  -> h[C_c] + amp[C_c] + so3[C_c] + trdox[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012420339_1');%original 'WP_012420339_1' aps[C_c] + trdrd[C_c]  -> 2 h[C_c] + amp[C_c] + so3[C_c] + trdox[C_c] 

model_AM = removeRxns(model_AM, 'R_BTS_1');
model_AM=addReaction(model_AM, 'R_BTS_1','reactionName', 'Biotin synthase', 'reactionformula', 'dtbt[C_c] + 2 s[C_c] 	->	h2s[C_c] + btn[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419545_1');%original dtbt[C_c] + 2 s[C_c] 	->	h[C_c] + h2s[C_c] + btn[C_c]

model_AM = removeRxns(model_AM, 'R_FMNRx');
model_AM=addReaction(model_AM, 'R_FMNRx','reactionName', 'FMN reductase', 'reactionformula', '2 h[C_c] + nadh[C_c] + fmn[C_c] -> nad[C_c] + fmnh2[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_042449258_1');%original h[C_c] + nadh[C_c] + fmn[C_c] 	->	nad[C_c] + fmnh2[C_c]

model_AM = removeRxns(model_AM, 'R_AKP1');
model_AM=addReaction(model_AM, 'R_AKP1','reactionName', 'Alkaline phosphatase  Dihydroneopterin', 'reactionformula', '3 h2o[C_c] + ahdt[C_c] -> 3 pi[C_c] + h[C_c] + dhnpt[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419918_1');%original 3 h2o[C_c] + ahdt[C_c] 	->	3 pi[C_c] + 2 h[C_c] + dhnpt[C_c]

model_AM = removeRxns(model_AM, 'R_FMNAT');
model_AM=addReaction(model_AM, 'R_FMNAT','reactionName', 'FMN adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + fmn[C_c] -> ppi[C_c] + fad[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419799_1');%original atp[C_c] + h[C_c] + fmn[C_c] 	->	ppi[C_c] + fad[C_c]

model_AM = removeRxns(model_AM, 'R_PRAIS');
model_AM=addReaction(model_AM, 'R_PRAIS','reactionName', 'Phosphoribosylaminoimidazole synthase', 'reactionformula', 'atp[C_c] + fpram[C_c] ->	adp[C_c] + pi[C_c] + h[C_c] + air[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419510_1');%,'geneRule','WP_012419510_1' original atp[C_c] + gly[C_c] + pram[C_c] 	<=>	adp[C_c] + pi[C_c] + h[C_c] + gar[C_c]

model_AM = removeRxns(model_AM, 'R_PANTS');
model_AM=addReaction(model_AM, 'R_PANTS','reactionName', 'Pantothenate synthase', 'reactionformula', 'atp[C_c] + ala_B[C_c] + pant__R[C_c] -> 2 h[C_c] + amp[C_c] + ppi[C_c] + pnto__R[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419457_1');%,'geneRule','WP_012419457_1' original atp[C_c] + ala_B[C_c] + pant__R[C_c] 	->	h[C_c] + amp[C_c] + ppi[C_c] + pnto__R[C_c]

model_AM = removeRxns(model_AM, 'R_PTPATi');
model_AM=addReaction(model_AM, 'R_PTPATi','reactionName', 'Pantetheine-phosphate adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + pan4p[C_c] 	->	ppi[C_c] + dpcoa[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012420770_1'); % ,'geneRule','WP_012420770_1' original atp[C_c] + h[C_c] + pan4p[C_c] 	->	ppi[C_c] + dpcoa[C_c]

model_AM = removeRxns(model_AM, 'R_PRAMPC_1');
model_AM=addReaction(model_AM, 'R_PRAMPC_1','reactionName', 'Phosphoribosyl AMP cyclohydrolase', 'reactionformula', 'h2o[C_c] + prbamp[C_c] -> prfp[C_c]','lowerBound',0,'upperBound',1000 ,'geneRule','WP_012419663_1');% ,'geneRule','WP_012419663_1' original h2o[C_c] + h[C_c] + prbamp[C_c]  -> prfp[C_c] 

model_AM = removeRxns(model_AM, 'R_PAPSR');
model_AM=addReaction(model_AM, 'R_PAPSR','reactionName', 'Phosphoadenylyl-sulfate reductase (thioredoxin)', 'reactionformula', 'paps[C_c] + trdrd[C_c]  -> h[C_c] + pap[C_c] + so3[C_c] + trdox[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419743_1');%,'geneRule','WP_012419743_1' original paps[C_c] + trdrd[C_c]  -> 2 h[C_c] + pap[C_c] + so3[C_c] + trdox[C_c]

model_AM = removeRxns(model_AM, 'R_PUNP4');
model_AM=addReaction(model_AM, 'R_PUNP4','reactionName', 'Purine-nucleoside phosphorylase (Deoxyguanosine)', 'reactionformula', 'pi[C_c] + dgsn[C_c] 	<=>	2dr1p[C_c] + gua[C_c] + h[C_c]','lowerBound',-1000,'upperBound',1000);%,'geneRule','WP_02419743_1' original pi[C_c] + dgsn[C_c] 	<=>	2dr1p[C_c] + gua[C_c]

model_AM = removeRxns(model_AM, 'R_NTD9');
model_AM=addReaction(model_AM, 'R_NTD9','reactionName', '5-nucleotidase (GMP)', 'reactionformula', 'h2o[C_c] + gmp[C_c] ->	h[C_c] + pi[C_c] + gsn[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419969_1');% ,'geneRule','WP_012419969_1' original h2o[C_c] + gmp[C_c] 	->	pi[C_c] + gsn[C_c]

model_AM = removeRxns(model_AM, 'R_GUAPRT');
model_AM=addReaction(model_AM, 'R_GUAPRT','reactionName', 'Guanine phosphoribosyltransferase', 'reactionformula', 'h[C_c] + prpp[C_c] + gua[C_c]  -> ppi[C_c] + gmp[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012420964_1');%,'geneRule','WP_012420964_1' original prpp[C_c] + gua[C_c]  -> ppi[C_c] + gmp[C_c]

model_AM = removeRxns(model_AM, 'R_IGPDH_1');
model_AM=addReaction(model_AM, 'R_IGPDH_1','reactionName', 'Imidazoleglycerol phosphate dehydratase', 'reactionformula', 'eig3p[C_c]  -> h2o[C_c] + imacp[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012420511_1');%,'geneRule','WP_012420511_1' original h[C_c] + eig3p[C_c]  -> h2o[C_c] + imacp[C_c] 

model_AM = removeRxns(model_AM, 'R_PRAGSr');
model_AM=addReaction(model_AM, 'R_PRAGSr','reactionName', 'Phosphoribosylglycinamide synthase', 'reactionformula', 'atp[C_c] + gly[C_c] + pram[C_c]  <=> adp[C_c] + pi[C_c] + 2 h[C_c] + gar[C_c]','lowerBound',-1000,'upperBound',1000,'geneRule','WP_012420890_1');%,'geneRule','WP_012420890_1' original atp[C_c] + gly[C_c] + pram[C_c]  <=> adp[C_c] + pi[C_c] + h[C_c] + gar[C_c] 

model_AM = removeRxns(model_AM, 'R_NADDP');
model_AM=addReaction(model_AM, 'R_NADDP','reactionName', 'NAD diphosphatase', 'reactionformula', 'h2o[C_c] + nad[C_c] 	->	3 h[C_c] + amp[C_c] + nmn[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419073_1');%,'geneRule','WP_012419073_1' original h2o[C_c] + nad[C_c] 	->	2 h[C_c] + amp[C_c] + nmn[C_c]

model_AM = removeRxns(model_AM, 'R_DNTPPA');
model_AM=addReaction(model_AM, 'R_DNTPPA','reactionName', '5-nucleotidase (GMP)', 'reactionformula', 'h2o[C_c] + ahdt[C_c] -> ppi[C_c] + dhpmp[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419136_1 or WP_012419160_1 or WP_012420743_1 or WP_052294492_1');%original h2o[C_c] + ahdt[C_c] 	->	h[C_c] + ppi[C_c] + dhpmp[C_c]

model_AM = removeRxns(model_AM, 'R_NMNAT');
model_AM=addReaction(model_AM, 'R_NMNAT','reactionName', 'Nicotinamide-nucleotide adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + nmn[C_c] -> nad[C_c] + ppi[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419136_1 or WP_012419444_1');% ,'geneRule','( WP_012419136_1 or WP_012419444_1 )' original atp[C_c] + h[C_c] + nmn[C_c] 	->	nad[C_c] + ppi[C_c]

model_AM = removeRxns(model_AM, 'R_HSTPTr');
model_AM=addReaction(model_AM, 'R_HSTPTr','reactionName', 'Histidinol phosphate transaminase', 'reactionformula', 'glu__L[C_c] + imacp[C_c] -> akg[C_c] + hisp[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012420512_1 or WP_012421161_1');% ,'geneRule','( WP_012420512_1 or WP_012421161_1 )' original glu__L[C_c] + imacp[C_c]  -> h[C_c] + akg[C_c] + hisp[C_c] 

model_AM = removeRxns(model_AM, 'R_GTPCI');
model_AM=addReaction(model_AM, 'R_GTPCI','reactionName', 'GTP cyclohydrolase I', 'reactionformula', 'h2o[C_c] + gtp[C_c]  -> 2 h[C_c] + for[C_c] + ahdt[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_012419153_1 or WP_012420914_1');%,'geneRule','( WP_012419153_1 or WP_012420914_1 )' original h2o[C_c] + gtp[C_c]  -> h[C_c] + for[C_c] + ahdt[C_c] 

model_AM = removeRxns(model_AM, 'R_OACT');
model_AM=addReaction(model_AM, 'R_OACT','reactionName', 'O-antigen Acetyl-Transferase', 'reactionformula', '2 h[C_c] + accoa[C_c] + udcpo4[C_c] -> coa[C_c] + udcpo5[C_c]','lowerBound',0,'upperBound',1000);%,'geneRule','' original accoa[C_c] + udcpo4[C_c] 	->	coa[C_c] + udcpo5[C_c]

model_AM = removeRxns(model_AM, 'R_SULR_syn');
model_AM=addReaction(model_AM, 'R_SULR_syn','reactionName', 'Ferredoxin-sulfite reductase', 'reactionformula', '7 h[C_c] + so3[C_c] + 6 fdxrd[C_c]  -> 3 h2o[C_c] + h2s[C_c] + 6 fdxo_2_2[C_c]','lowerBound',0,'upperBound',1000);%,'geneRule','' original 8 h[C_c] + so3[C_c] + 6 fdxrd[C_c]  -> 3 h2o[C_c] + h2s[C_c] + 6 fdxo_2_2[C_c] 


% Remove duplicated reactions and loop reactions
model_AM = removeRxns(model_AM,{'R_ACONTa','R_ACONTb','R_CO2tex','R_CO2tpp','R_H2Ot','R_NH4tex','R_NH4tpp','R_NO3t7pp','R_GLYAT','R_ICH','R_NTP1','R_PGK_1','R_PGM_1','R_PRAIS','R_TDPDRR','R_BTS_1','R_CPPPGOAN2_1','R_IGPDH_1','R_PRAMPC_1','R_KARA1','R_DHFOR2','R_FOMETRi','R_METS_1','R_SULR'}); % NTP1 was atp[C_c] + h2o[C_c] => adp[C_c] + pi[C_c] + h[C_c];check with MetaCyc 

% Update reactions bounds
model_AM = changeRxnBounds(model_AM,{'R_GK1','R_GK2','R_ADK3','R_ADK4','R_NDPK1','R_NDPK8','R_DADK','R_NDPK9', 'R_FRD2rpp'},[0,0,0,0,0,0,0,0,-1000],'l');
model_AM = changeRxnBounds(model_AM, 'R_FRD2rpp',0,'u');% 
model_AM = changeRxnBounds(model_AM, 'R_FERIRDe',0,'b');% 
model_AM = changeRxnBounds(model_AM,'R_PPAt4pp',-1000,'l');% From Metacyc https://biocyc.org/FAECPRAU/class-tree?object=Transport-Reactions

% Reactions modified because unlikely thermodynamically feasible
    
    model_AM = changeRxnBounds(model_AM,{'R_NDPK1';'R_NDPK2';'R_NDPK3';'R_NDPK4';'R_NDPK5';'R_NDPK6';'R_NDPK7';'R_NDPK8'},repelem(0,8),'l'); %Nucleoside-diphosphatekinase
    model_AM = changeRxnBounds(model_AM,{'R_HSDx';'R_HSDy'},repelem(0,2),'u');%L-homoserine dehydrogenase https://metacyc.org/gene?orgid=META&id=ASPKINIIHOMOSERDEHYDROGII-MONOMER
 

%% Check for duplicated reactions (method 'FR')
checkDuplicateRxn(model_AM,'FR',1,1);% no duplicates found

%% Change the ATPM name for 
model_AM = removeRxns(model_AM, {'R_ATPM'});
model_AM=addReaction(model_AM, 'R_ngam','reactionName', 'Non-growth_ATP_maintenance', 'reactionformula', 'atp[C_c] + h2o[C_c]  -> adp[C_c] + pi[C_c] + h[C_c]','geneRule','spontaneous','lowerBound',0,'upperBound',1000);

%% Reactions added to the model
    
        % Add metabolites 
        
        Met={'acgal[C_e]','acgal[C_c]','acgal6p[C_c]','acgam[C_e]','asp__L[C_e]','citr__L[C_e]','galam6p[C_c]','gam[C_e]','gln__L[C_e]','his__L[C_e]','lys__L[C_e]','mal__L[C_e]','mal__L[C_p]','met__L[C_e]','phe__L[C_e]','ser__L[C_e]','succ[C_e]','trp__L[C_e]','tyr__L[C_e]','fol[C_e]'};

    model_AM= addMultipleMetabolites(model_AM,Met);
    
        % Add Exchange reactions
        rxnToAdd.metList={'acgal[C_e]','acgam[C_e]','asp__L[C_e]','citr__L[C_e]','gam[C_e]','gln__L[C_e]','his__L[C_e]','mal__L[C_e]','met__L[C_e]','lys__L[C_e]','phe__L[C_e]','ser__L[C_e]','na1[C_e]','succ[C_e]','trp__L[C_e]','tyr__L[C_e]','fol[C_e]'};
        rxnToAdd.lb    = [repelem(-1000,length(rxnToAdd.metList))];
        rxnToAdd.ub    = [repelem(1000,length(rxnToAdd.metList))];
        
    model_AM = addExchangeRxn_v2(model_AM,rxnToAdd.metList,rxnToAdd.lb,rxnToAdd.ub); %the function has been modified to add 'R_' in front of exchange reaction names

        % Add transport reactions 
        
    model_AM=addReaction(model_AM, 'R_SUCCt2r','reactionName', 'Succinate transport via proton symport', 'reactionformula', 'h[C_e] + succ[C_e] <=> h[C_c] + succ[C_c]'); %Thiele I, et al (2017) model
    model_AM=addReaction(model_AM, 'R_ACGApts','reactionName', 'N-Acetyl-D-glucosamine transport via PEP:Pyr PTS', 'reactionformula', 'acgam[C_e] + pep[C_c] -> acgam6p[C_c] + pyr[C_c] '); %Thiele I, et al (2017) model
    model_AM=addReaction(model_AM, 'R_ACGALt','reactionName', 'N-Acetyl-D-galactosamine diffusion', 'reactionformula', 'acgal[C_e] + h[C_e] <=> acgal[C_c] + h[C_c]'); %Thiele I, et al (2017) model
    model_AM=addReaction(model_AM, 'R_GAMpts','reactionName', 'D-glucosamine transport via PEP:Pyr PTS', 'reactionformula', 'pep[C_c] + gam[C_e] -> pyr[C_c] + gam6p[C_c]'); %http://bigg.ucsd.edu/universal/reactions/GAMpts
    model_AM=addReaction(model_AM, 'R_MALt5','reactionName', 'Na+/malate symporter', 'reactionformula', 'mal__L[C_e] + na1[C_e] <=> mal__L[C_p] + na1[C_p]'); %http://bigg.ucsd.edu/universal/reactions/MALt5
    model_AM=addReaction(model_AM, 'R_MALtp5','reactionName', 'Malate diffusion', 'reactionformula', 'mal__L[C_p] <=> mal__L[C_c]'); %assumption
    model_AM=addReaction(model_AM, 'R_ACGALK3','reactionName', 'N-Acetyl-D-galactosamine kinase', 'reactionformula', 'atp[C_c] + acgal[C_c] -> adp[C_c] + acgal6p[C_c] + h[C_c]'); %
    model_AM=addReaction(model_AM, 'R_ACGAL6P_deacetylase','reactionName', 'N-Acetyl-D-galactosamine 6 phosphate deacetylase', 'reactionformula', 'acgal6p[C_c] + h2o[C_c] + h[C_c] -> ac[C_c] + galam6p[C_c]'); %
    model_AM=addReaction(model_AM, 'R_GALAM6P_deacetylase','reactionName', 'Galactosamine-6-phosphate deaminase', 'reactionformula', 'galam6p[C_c] + h2o[C_c] -> nh4[C_c] + tag6p__D[C_c] + h[C_c]'); %

            % Amino Acids transpoters (Assumed symporter)
   
   model_AM=addReaction(model_AM, 'R_ASPt2r','reactionName', 'Aspartate reversible transport via proton symport', 'reactionformula', 'h[C_e] + asp__L[C_e] <=> h[C_c] + asp__L[C_c]'); 
   model_AM=addReaction(model_AM, 'R_Cystr','reactionName', 'L-Cysteine reversible transport via proton symport', 'reactionformula', 'h[C_e] + cys__L[C_e] <=> h[C_c] + cys__L[C_c]'); 
   model_AM=addReaction(model_AM, 'R_GLNt2r','reactionName', 'Glutamine reversible transport via proton symport', 'reactionformula', 'h[C_e] + gln__L[C_e] <=> h[C_c] + gln__L[C_c]');  
   model_AM=addReaction(model_AM, 'R_GLUt2r','reactionName', 'Glutamate reversible transport via proton symport', 'reactionformula', 'h[C_e] + glu__L[C_e] <=> h[C_c] + glu__L[C_c]'); 
   model_AM=addReaction(model_AM, 'R_Hist2r','reactionName', 'L-Histidine reversible transport via proton symport', 'reactionformula', 'h[C_e] + his__L[C_e] <=> h[C_c] + his__L[C_c]');
   model_AM=addReaction(model_AM, 'R_Lyst2r','reactionName', 'L-Lysine reversible transport via proton symport', 'reactionformula', 'h[C_e] + lys__L[C_e] <=> h[C_c] + lys__L[C_c]');
   model_AM=addReaction(model_AM, 'R_Met2r','reactionName', 'L-Methionine reversible transport via proton symport', 'reactionformula', 'h[C_e] + met__L[C_e] <=> h[C_c] + met__L[C_c]');
   model_AM=addReaction(model_AM, 'R_PHEt2r','reactionName', 'L-Phenylalanine reversible transport via proton symport', 'reactionformula', 'h[C_e] + phe__L[C_e] <=> h[C_c] + phe__L[C_c]'); 
   model_AM=addReaction(model_AM, 'R_SERt2r','reactionName', 'L-Serine reversible transport via proton symport', 'reactionformula', 'h[C_e] + ser__L[C_e] <=> h[C_c] + ser__L[C_c]'); 
   model_AM=addReaction(model_AM, 'R_TRPt2r','reactionName', 'L-Tryptophan reversible transport via proton symport', 'reactionformula', 'h[C_e] + trp__L[C_e] <=> h[C_c] + trp__L[C_c]');
   model_AM=addReaction(model_AM, 'R_CITRt2r','reactionName', 'L-Citrulline reversible transport via proton symport', 'reactionformula', 'h[C_e] + citr__L[C_e] <=> h[C_c] + citr__L[C_c]');
   model_AM=addReaction(model_AM, 'R_Tyrt2r','reactionName', 'L-Tyrosine reversible transport via proton symport', 'reactionformula', 'h[C_e] + tyr__L[C_e] <=> h[C_c] + tyr__L[C_c]');
            

        % Folate transporter 
      model_AM=addReaction(model_AM, 'R_Foltr','reactionName', 'Folate transport via proton symport', 'reactionformula', 'h[C_e] + fol[C_e] -> h[C_c] + fol[C_c]'); %BT_3620 https://www.uniprot.org/uniprot/Q8A1N9
      model_AM=addReaction(model_AM, 'R_FRD4rpp','reactionName', 'Fumarate reductase / succinate dehydrogenase (irreversible) (periplasmic, membrane potential dissipating)', 'reactionformula', ' fum[C_c] + qh2[C_c] -> succ[C_c] + q[C_c]','lowerBound',0,'upperBound',1000);%,'geneRule','' original 8 h[C_c] + so3[C_c] + 6 fdxrd[C_c]  -> 3 h2o[C_c] + h2s[C_c] + 6 fdxo_2_2[C_c] 

 

%% Add formulas and Names of metabolites added before
    
        Names={'N-Acetyl-D-galactosamine','N-Acetyl-D-galactosamine','N-Acetyl-D-galactosamine-6-phosphate','N-Acetyl-D-glucosamine','L-Aspartate','L-Citrulline','D-Galactosamine-6-phosphate','D-Glucosamine','L-Glutamine','L-Histidine','L-Lysine','L-Malate','L-Malate','L-Methionine','L-Phenylalanine','L-Serine','Succinate','L-Tryptophan','L-Tyrosine','Folate'};
        Formulas={'C8H15NO6','C8H15NO6','C8H14NO9P','C8H15NO6','C4H6NO4','C6H13N3O3','C6H14NO8P','C6H14NO5','C5H10N2O3','C6H9N3O2','C6H15N2O2','C4H4O5','C4H4O5','C5H11NO2S','C9H11NO2','C3H7NO3','C4H4O4','C11H12N2O2','C9H11NO3','C19H18N7O6'};
    
        for i=1:numel(Met)            
            pos = strcmp(model_AM.mets,char(Met(i)));
            k= find(pos);
            model_AM.metNames(k)= Names(i);
            model_AM.metFormulas(k) = Formulas(i);
        end


 %% Model constrains
    
            %Amino Acids
            AA= {'R_EX_ala__L_e';'R_EX_arg__L_e';'R_EX_asn__L_e';'R_EX_asp__L_e';'R_EX_citr__L_e';'R_EX_cys__L_e';'R_EX_gln__L_e';'R_EX_glu__L_e';'R_EX_gly_e';'R_EX_his__L_e';'R_EX_ile__L_e';'R_EX_leu__L_e';'R_EX_lys__L_e';'R_EX_met__L_e';'R_EX_phe__L_e';'R_EX_pro__L_e';'R_EX_ser__L_e';'R_EX_thr__L_e';'R_EX_trp__L_e';'R_EX_tyr__L_e';'R_EX_val__L_e'};
            model_AM = changeRxnBounds(model_AM,AA,[repelem(0,length(AA))],'u'); 
            model_AM = changeRxnBounds(model_AM,AA,[repelem(-1000,length(AA))],'l'); 
    
            AA_other={'R_EX_ala__D_e';'R_EX_alaala_e';'R_EX_ser__D_e'};
            model_AM = changeRxnBounds(model_AM,AA_other,[repelem(0,length(AA_other))],'u'); 
            model_AM = changeRxnBounds(model_AM,AA_other,[repelem(0,length(AA_other))],'l'); 
           
            % Nucleic Acids
          
            NA_other= {'R_EX_gua_e';'R_EX_thym_e';'R_EX_cmp_e'};
            model_AM = changeRxnBounds(model_AM,NA_other,[repelem(0,length(NA_other))],'u'); 
            model_AM = changeRxnBounds(model_AM,NA_other,[repelem(0,length(NA_other))],'l');
            
            % Carbohydrates
            CH= {'R_EX_abt__L_e';'R_EX_anhgm_e';'R_EX_arab__D_e';'R_EX_arab__L_e';'R_EX_dha_e';'R_EX_fru_e';'R_EX_galt_e';'R_EX_galur_e';'R_EX_gam6p_e';'R_EX_glc__D_e';'R_EX_glyc_e';'R_EX_lcts_e';'R_EX_man_e';'R_EX_melib_e';'R_EX_raffin_e';'R_EX_rbt_e';'R_EX_sbt__D_e'};
            model_AM = changeRxnBounds(model_AM,CH,[repelem(0,length(CH))],'b');

            % Mucin Component (carbohydrate part)
            MCP= {'R_EX_acgal_e';'R_EX_acgam_e';'R_EX_acnam_e';'R_EX_fuc__L_e';'R_EX_gal_e';'R_EX_gam_e'}; %According to Thiele mannose is not in mucin 
            model_AM = changeRxnBounds(model_AM,MCP,[-1000,-1000,0,-1000,-1000,0],'l');
            model_AM = changeRxnBounds(model_AM,MCP,[repelem(0,6)],'u');

            %Ions and Vitamins
            IV={'R_EX_ca2_e';'R_EX_cl_e';'R_EX_cobalt2_e';'R_EX_cu2_e';'R_EX_fe2_e';'R_EX_fe3_e';'R_EX_h_e';'R_EX_k_e';'R_EX_mg2_e';'R_EX_mn2_e';'R_EX_na1_e';'R_EX_nh4_e';'R_EX_no2_e';'R_EX_pi_e';'R_EX_pnto__R_e';'R_EX_so4_e';'R_EX_zn2_e'};% cobalt2, cu2, fe2, fe3, zn2, panthotenate (pnto), pyridoxal (pydx) and riboflavin (ribfl) are not added in the medium but can be in Yeast Extract
            model_AM = changeRxnBounds(model_AM,IV,[repelem(1000,length(IV))],'u'); 
            model_AM = changeRxnBounds(model_AM,IV,[repelem(-1000,length(IV))],'l'); 

            % Gaz and water
            Gaz={'R_EX_co2_e';'R_EX_h2_e';'R_EX_h2o_e';'R_EX_h2s_e';'R_EX_o2_e'};
            model_AM = changeRxnBounds(model_AM,Gaz,[1000,0,1000,0,0],'u'); 
            model_AM = changeRxnBounds(model_AM,Gaz,[0,0,-1000,0,0],'l'); 
          
            % Lipid and Membrane components
            LipMc={'R_EX_etha_e';'R_EX_glcur_e';'R_EX_gm1lipa_e';'R_EX_gm2lipa_e';'R_EX_udcpo5_e';'R_EX_kdo2lipid4_e'};
            model_AM = changeRxnBounds(model_AM,LipMc,[repelem(0,length(LipMc))],'b');
                        
            % Siderophore (enterobactin)
            Sid={'R_EX_fe3hox_e';'R_EX_fe3hox_un_e';'R_EX_feenter_e';'R_EX_salchs4_e';'R_EX_pheme_e'};
            model_AM = changeRxnBounds(model_AM,Sid,[repelem(0,length(Sid))],'b');

            % Required fluxes to optimize Biomass
            model_AM = changeRxnBounds(model_AM,{'R_EX_glyclt_e'},1,'u'); 
            model_AM = changeRxnBounds(model_AM,{'R_EX_glyclt_e'},0,'l'); %if bz=0 or/and tton=0 => No growth;
       
            % Unclassified
            model_AM = changeRxnBounds(model_AM,'R_EX_acald_e',0,'u'); 
            model_AM = changeRxnBounds(model_AM,'R_EX_acald_e',-1000,'l'); 
            model_AM = changeRxnBounds(model_AM,'R_EX_no2_e',0,'l');
            model_AM = changeRxnBounds(model_AM,'R_EX_no3_e',0,'l'); 

            Unclass={'R_EX_14glucan_e';'R_EX_2pglyc_e';'R_EX_4abut_e';'R_EX_4hthr_e';'R_EX_4hphac_e';'R_EX_butso3_e';'R_EX_chol_e';'R_EX_chols_e';'R_EX_chtbs_e';'R_EX_ethso3_e';'R_EX_hxan_e';'R_EX_id3acald_e';'R_EX_indole_e';'R_EX_inost_e';'R_EX_isetac_e';'R_EX_lmn2_e';'R_EX_lmn30_e';'R_EX_meoh_e';'R_EX_mso3_e';'R_EX_pyr_e';'R_EX_s_e';'R_EX_spmd_e';'R_EX_sulfac_e';'R_EX_tartr__D_e';'R_EX_tma_e';'R_EX_tmao_e';'R_EX_tol_e'};
            model_AM = changeRxnBounds(model_AM,Unclass,[repelem(0,length(Unclass))],'b');
                
            % Experimental rates
                
            model_AM = changeRxnBounds(model_AM, 'R_EX_ac_e',10.8,'u');
            model_AM = changeRxnBounds(model_AM, 'R_EX_ac_e',9.21,'l');
            
            model_AM = changeRxnBounds(model_AM, 'R_EX_cit_e',0.79,'u');
            model_AM = changeRxnBounds(model_AM, 'R_EX_cit_e',-0.38,'l');

            model_AM = changeRxnBounds(model_AM,'R_EX_succ_e',0,'u');
            model_AM =changeRxnBounds(model_AM,'R_EX_succ_e',-0.3,'l');

            model_AM = changeRxnBounds(model_AM, 'R_EX_ppa_e',8.58,'u');
            model_AM = changeRxnBounds(model_AM, 'R_EX_ppa_e',7.38,'l');

            model_AM = changeRxnBounds(model_AM, 'R_EX_lac__D_e',0.134,'u'); 
            model_AM = changeRxnBounds(model_AM, 'R_EX_lac__D_e',0,'l'); 
          
            model_AM = changeRxnBounds(model_AM, 'R_EX_for_e',1.15,'u');
            model_AM = changeRxnBounds(model_AM, 'R_EX_for_e',0.59,'l');                  

            model_AM = changeRxnBounds(model_AM, 'R_EX_etoh_e',0,'u');
            model_AM = changeRxnBounds(model_AM, 'R_EX_etoh_e',0,'l');
            
            model_AM = changeRxnBounds(model_AM, 'R_EX_mal__L_e',-0.32,'u'); 
            model_AM = changeRxnBounds(model_AM, 'R_EX_mal__L_e',-0.75,'l'); 

            model_AM = changeRxnBounds(model_AM, 'R_EX_but_e',0,'u');
            model_AM = changeRxnBounds(model_AM, 'R_EX_but_e',0,'l');
                
 %% Biomass reaction modification.
     % Biomass GAM modified to reach experimental growth rate
     
     model_AM = addReaction(model_AM,'Growth','reactionName', 'Biomass reaction','reactionFormula','0.000215957 10fthf[C_c] + 26 atp[C_c] + 0.1 amp[C_c] + 0.254849 glu__L[C_c] + 26 h2o[C_c] + 0.00177318 nad[C_c] + 0.000557809 coa[C_c] + 0.000432883 nadp[C_c] + 0.124 gmp[C_c] + 0.000215957 amet[C_c] + 0.436303 leu__L[C_c] + 0.233442 asp__L[C_c] + 0.000215957 fad[C_c] + 0.497466 ala__L[C_c] + 0.00504062 cl[C_c] + 0.000215957 ribflv[C_c] + 0.254849 gln__L[C_c] + 0.0886878 cys__L[C_c] + 0.208977 ser__L[C_c] + 0.245675 thr__L[C_c] + 0.0917461 his__L[C_c] + 0.000215957 thmpp[C_c] + 0.074 cmp[C_c] + 0.0136 damp[C_c] + 0.0143 dgmp[C_c] + 0.59329 gly[C_c] + 0.148832 met__L[C_c] + 0.000215957 thf[C_c] + 0.233442 asn__L[C_c] + 0.000215957 pydx5p[C_c] + 0.286451 arg__L[C_c] + 9.68419e-05 mql8[C_c] + 0.214074 pro__L[C_c] + 0.332324 lys__L[C_c] + 0.00650293 fe2[C_c] + 0.179414 phe__L[C_c] + 0.004201 so4[C_c] + 0.189029 k[C_c] + 0.281354 ile__L[C_c] + 0.000686609 cu2[C_c] + 0.00504062 ca2[C_c] + 0.00840103 mg2[C_c] + 0.000669178 mn2[C_c] + 9.68419e-05 cobalt2[C_c] + 0.000330231 zn2[C_c] + 0.00756142 fe3[C_c] + 0.0802 ump[C_c] + 0.0137 dcmp[C_c] + 0.000215957 mlthf[C_c] + 0.0135 dtmp[C_c] + 0.0550478 trp__L[C_c] + 0.409798 val__L[C_c] + 0.133541 tyr__L[C_c] + 0.0968419 uaagmda[C_c] 	->	26 adp[C_c] + 26 pi[C_c] + 26 h[C_c]')
       
%% Optimization for biomass

    model_AM=changeObjective(model_AM,'Growth'); 
    FBAsolutionAM=optimizeCbModel(model_AM,'max','one',false); 
    printFluxVector(model_AM, FBAsolutionAM.x, 'true');

%% Export model as SMBL containing constrains

    writeCbModel(model_AM,'sbml','[your_path]/Akkermansia_muciniphila_ATCC_BAA_835',{'C_c','C_p','C_e'},{'cytosol','periplasm','extracellular space'}) %SBML level = 2 by default, sbmlVersion = 1 by default

%% Random Sampling (RAVEN)

        cd '[your_path]';
        model_AMrav=importModel('Akkermansia_muciniphila_ATCC_BAA_835.xml') 

    % Fix measured exchange fluxes around 5% of the value from FBA
    exIdx = getIndexes(model_AMrav,{'EX_ac_e','EX_cit_e','EX_succ_e','EX_ppa_e','EX_lac__D_e','EX_for_e','EX_etoh_e','EX_mal__L_e', 'EX_but_e','Growth','ngam'},'rxns');%

    fluxes = sol.x(exIdx);
    model_AMrav = setParam(model_AMrav,'var',exIdx,fluxes,10);    
    [~, goodrxn] = randomSampling(model_AMrav,1,true,true,true); % 5000 Sampling generated
    rs = randomSampling(model_AMrav,5000,true,true,true,goodrxn); % 5000 Sampling generated
    fluxMean = full(mean(rs,2));
    fluxSD = full(std(rs,0,2));
    
    %% Save FBA predicted fluxes
    
    cd '[your_path]';
 
    out=table(model_AMrav.rxns,model_AMrav.rxnNames,constructEquations(model_AMrav),sol.x);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Flux'};
    out(out.Flux==0,:) =[];
    writetable(out,'AM_FBA.xlsx');


    %% Save FVA predicted fluxes (mean+SD) and random sampling matrix
    clear out;out = fluxMean;
    clear out2;out2 = fluxSD;
    
    cd '[your_path]';

    out =table(model_AMrav.rxns,model_AMrav.rxnNames, constructEquations(model_AMrav), out ,out2);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Mean_RS' 'SD_RS'};
    out(out.Mean_RS==0,:) =[];
    writetable(out,'AM_RS.xlsx');

    rs(~any(rs,2),:) =[];
    RS_df=   array2table(full(rs));
    RS_df = [out(:,[1:3]),RS_df];
    writetable(RS_df,'AM_RS_5000samples.xlsx');

%% NADH and NADPH reactions
i= 'Nicotinamide adenine dinucleotide phosphate - reduced';
    [fluxes, rxnIdx] = getMetProduction(model_AMrav,i,fluxMean,true);   
    clear out;
    out.rxns    =model_AMrav.rxns(rxnIdx);
    out.rxnNames= model_AMrav.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model_AMrav,rxnIdx);
    out.fluxes  = num2cell(fluxes); 
    out = table(out.rxns,out.rxnNames,out.rxnEqns,out.fluxes);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Flux'};
    writetable(out,'AM_NADPH.xlsx');
    
i= 'Nicotinamide adenine dinucleotide - reduced';
    [fluxes, rxnIdx] = getMetProduction(model_AMrav,i,fluxMean,true);
    clear out;
    out.rxns    =model_AMrav.rxns(rxnIdx);
    out.rxnNames= model_AMrav.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model_AMrav,rxnIdx);
    out.fluxes  = num2cell(fluxes); 
    out = table(out.rxns,out.rxnNames,out.rxnEqns,out.fluxes);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Flux'};
    writetable(out,'AM_NADH.xlsx');

    
     %% Individual depletion analysis
    EX_RS= out(contains(out.(1),'EX_'),:);%48 rows
    temp= EX_RS(EX_RS.Mean_RS <0,:);
    IV2=strrep(IV,'R_','');
    idx = ismember(temp.Reaction,IV2); %get indexes of ions and vitamins metabolites
    temp(idx,:)=[];%27
    
    % Remove exp metabolites
    Exp_meta={'EX_ac_e','EX_but_e','EX_cit_e','EX_succ_e','EX_for_e','EX_etoh_e','EX_ppa_e','EX_lac__D_e','EX_mal__L_e'};
    idx2 = ismember(temp.Reaction,Exp_meta); %get indexes of ions and vitamins metabolites
    temp(idx2,:)=[];%26
    
    % Remove water
    temp(strcmp(temp.Reaction,'EX_h2o_e'),:)=[]; %25
    
    temp.Reaction = strcat('R_', temp.Reaction);
    temp.name = strrep(temp.Reaction,'R_EX_','');
    
    cd '[your_path]';
    GR_prd = table();
    
    for i=1:length(temp.Reaction)        
        model_AM = changeRxnBounds(model_AM, temp.Reaction(i),0,'b');% 
        FBAsolution_temp=optimizeCbModel(model_AM,'max','one',false); % Optimize for GR
        disp([temp.Reaction(i) FBAsolution_temp.f])
        GR_prd.GR(i)= FBAsolution_temp.f;
        GR_prd.Depelted_Metabolite(i)= temp.Reaction(i);
            if isnan(FBAsolution_temp.f)~=1
                out=table(model_AM.rxns,model_AM.rxnNames,FBAsolution_temp.x);
                out.Properties.VariableNames = {'Reaction' 'Description' 'Flux'};
                out(out.Flux==0,:) =[];
                writetable(out,string(strcat('AM_FBA_WO_',temp.Reaction(i),'.xlsx')));
            end
        model_AM = changeRxnBounds(model_AM,  temp.Reaction(i),-1000,'l');%
    end
writetable(GR_prd,string(strcat('AM_FBA_GR_indDepletion.xlsx')));
    