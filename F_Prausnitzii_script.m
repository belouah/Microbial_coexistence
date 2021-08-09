%initialization CobraToolbox and setting GUROBI as solver
initCobraToolbox
changeCobraSolver('gurobi','all')

%% Load the model

getFullPath(fullfile(cd, 'Faecalibacterium_prausnitzii_A2_165.xml'));
cd '[your_path]';
model=readCbModel('Faecalibacterium_prausnitzii_A2_165','fileType','SBML');

%% Rename the model
model_FP = model;

%% Curation
%Update metabolites formula

pos1 = strcmp(model_FP.metNames,'L-Cystathionine');
k= find(pos1);
model_FP.metFormulas(k);
model_FP.metFormulas(k)= {'C7H14N2O4S'};

pos2 = strcmp(model_FP.metNames,'Pyridoxal 5''-phosphate');
k2= find(pos2);
model_FP.metFormulas(k2);
model_FP.metFormulas(k2)= {'C8H8NO6P'};

pos4 = strcmp(model_FP.metNames,'2-Amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl)dihydropteridine triphosphate');
k4= find(pos4);
model_FP.metFormulas(k4);
model_FP.metFormulas(k4)= {'C9H12N5O13P3'};

% Update existing reactions
model_FP = removeRxns(model_FP, {'R_FEENTERR2'});
model_FP=addReaction(model_FP, 'R_FEENTERR2','reactionName', 'Fe-enterobactin reduction (Fe(III)-unloading)', 'reactionformula', 'fmnh2[C_c] + 2 feenter[C_c] -> 3 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 enter[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original fmnh2[C_c] + 2 feenter[C_c] 	->	2 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 enter[C_c]

model_FP = removeRxns(model_FP, {'R_FMNAT'});
model_FP=addReaction(model_FP, 'R_FMNAT','reactionName', 'FMN adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + fmn[C_c]  -> ppi[C_c] + fad[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_005933136_1');%original atp[C_c] + h[C_c] + fmn[C_c]  -> ppi[C_c] + fad[C_c] 

model_FP = removeRxns(model_FP, {'R_FMNRx'});
model_FP=addReaction(model_FP, 'R_FMNRx','reactionName', 'FMN reductase', 'reactionformula', '2 h[C_c] + nadh[C_c] + fmn[C_c] -> nad[C_c] + fmnh2[C_c]','geneRule','WP_005931883_1','lowerBound',0,'upperBound',1000);%original h[C_c] + nadh[C_c] + fmn[C_c] 	->	nad[C_c] + fmnh2[C_c]

model_FP = removeRxns(model_FP, {'R_FMNRy'});
model_FP=addReaction(model_FP, 'R_FMNRy','reactionName', 'FMN reductase  NADPH dependent', 'reactionformula', '2 h[C_c] + nadph[C_c] + fmn[C_c] -> nadp[C_c] + fmnRD[C_c]','geneRule','WP_005935658_1','lowerBound',0,'upperBound',1000);%original h[C_c] + nadph[C_c] + fmn[C_c] 	->	nadp[C_c] + fmnRD[C_c]

model_FP = removeRxns(model_FP, {'R_PRAGSr'});
model_FP=addReaction(model_FP, 'R_PRAGSr','reactionName', 'Phosphoribosylglycinamide synthase', 'reactionformula', 'atp[C_c] + gly[C_c] + pram[C_c] <=>	adp[C_c] + pi[C_c] + 2 h[C_c] + gar[C_c]','geneRule','WP_005932200_1','lowerBound',-1000,'upperBound',1000);%original atp[C_c] + gly[C_c] + pram[C_c] 	<=>	adp[C_c] + pi[C_c] + h[C_c] + gar[C_c]

model_FP = removeRxns(model_FP, {'R_MTHFR3'});
model_FP=addReaction(model_FP, 'R_MTHFR3','reactionName', '5 10 methylenetetrahydrofolatereductase  NADPH', 'reactionformula', 'h[C_c] + nadph[C_c] + mlthf[C_c] ->	nadp[C_c] + 5mthf[C_c]','geneRule','WP_005936540_1','lowerBound',0,'upperBound',1000);%original 2 h[C_c] + nadph[C_c] + mlthf[C_c] 	->	nadp[C_c] + 5mthf[C_c]

model_FP = removeRxns(model_FP, {'R_PUNP4'});
model_FP=addReaction(model_FP, 'R_PUNP4','reactionName', 'Purine-nucleoside phosphorylase (Deoxyguanosine)', 'reactionformula', 'pi[C_c] + dgsn[C_c] 	<=>	2dr1p[C_c] + gua[C_c] + h[C_c]','geneRule','WP_035393670_1','lowerBound',-1000,'upperBound',1000);%original pi[C_c] + dgsn[C_c] 	<=>	2dr1p[C_c] + gua[C_c]

model_FP = removeRxns(model_FP, {'R_PTPATi'});
model_FP=addReaction(model_FP, 'R_PTPATi','reactionName', 'Pantetheine-phosphate adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + pan4p[C_c] 	->	ppi[C_c] + dpcoa[C_c]','geneRule','WP_005932740_1','lowerBound',0,'upperBound',1000);%original atp[C_c] + h[C_c] + pan4p[C_c] 	->	ppi[C_c] + dpcoa[C_c]

model_FP = removeRxns(model_FP, {'R_GUAPRT'});
model_FP=addReaction(model_FP, 'R_GUAPRT','reactionName', 'Guanine phosphoribosyltransferase', 'reactionformula', 'h[C_c] + prpp[C_c] + gua[C_c] ->	ppi[C_c] + gmp[C_c]','geneRule','( WP_005931845_1 or WP_005933386_1 )','lowerBound',0,'upperBound',1000);%original prpp[C_c] + gua[C_c] 	->	ppi[C_c] + gmp[C_c]

model_FP = removeRxns(model_FP, {'R_UDCPAT'});
model_FP=addReaction(model_FP, 'R_UDCPAT','reactionName', 'Abequosyl transferase', 'reactionformula', 'cdpabeq[C_c] + udcpgrm[C_c] -> 3 h[C_c] + cdp[C_c] + udcpo4[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original cdpabeq[C_c] + udcpgrm[C_c] 	->	h[C_c] + cdp[C_c] + udcpo4[C_c]

model_FP = removeRxns(model_FP, {'R_OACT'});
model_FP=addReaction(model_FP, 'R_OACT','reactionName', 'O-antigen Acetyl-Transferase', 'reactionformula', '2 h[C_c] + accoa[C_c] + udcpo4[C_c] -> coa[C_c] + udcpo5[C_c]','geneRule','','lowerBound',0,'upperBound',1000);%original accoa[C_c] + udcpo4[C_c] 	->	coa[C_c] + udcpo5[C_c]

model_FP = removeRxns(model_FP, {'R_NTD9'});
model_FP=addReaction(model_FP, 'R_NTD9','reactionName', '5''-nucleotidase (GMP)', 'reactionformula', 'h2o[C_c] + gmp[C_c]  -> pi[C_c] + gsn[C_c] + h[C_c]','lowerBound',0,'upperBound',1000,'geneRule','( WP_005932613_1 or WP_005935454_1 )');%original h2o[C_c] + gmp[C_c]  -> pi[C_c] + gsn[C_c] 

model_FP = removeRxns(model_FP, {'R_PNP'});
model_FP=addReaction(model_FP, 'R_PNP','reactionName', 'Histidinol phosphate transaminase', 'reactionformula', 'pi[C_c] + rnam[C_c]  <=> ncam[C_c] + r1p[C_c]','lowerBound',-1000,'upperBound',1000,'geneRule','WP_035393670_1');%original pi[C_c] + rnam[C_c]  <=> h[C_c] + ncam[C_c] + r1p[C_c] 

model_FP = removeRxns(model_FP, {'R_SALCHS4FER2'});
model_FP=addReaction(model_FP, 'R_SALCHS4FER2','reactionName', 'Salmochelin S4 Fe III reduction Fe III unloading', 'reactionformula', 'fmnh2[C_c] + 2 salchs4fe[C_c]  -> 3 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 salchs4[C_c]','lowerBound',0,'upperBound',1000,'geneRule','');%original fmnh2[C_c] + 2 salchs4fe[C_c]  -> 2 h[C_c] + fmn[C_c] + 2 fe2[C_c] + 2 salchs4[C_c]

model_FP = removeRxns(model_FP, {'R_NMNAT'});
model_FP=addReaction(model_FP, 'R_NMNAT','reactionName', 'Nicotinamide-nucleotide adenylyltransferase', 'reactionformula', 'atp[C_c] + 2 h[C_c] + nmn[C_c]  -> nad[C_c] + ppi[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_005932886_1');%original atp[C_c] + h[C_c] + nmn[C_c]  -> nad[C_c] + ppi[C_c] 

model_FP = removeRxns(model_FP, {'R_NMNDA'});
model_FP=addReaction(model_FP, 'R_NMNDA','reactionName', 'Nicotinamide-nucleotide amidase', 'reactionformula', 'h[C_c] + h2o[C_c] + nmn[C_c]  -> nh4[C_c] + nicrnt[C_c]','lowerBound',0,'upperBound',1000,'geneRule','');%original h2o[C_c] + nmn[C_c]  -> nh4[C_c] + nicrnt[C_c] 

model_FP = removeRxns(model_FP, {'R_TTONR1pp'});
model_FP=addReaction(model_FP, 'R_TTONR1pp','reactionName', 'Trithionate Reductase', 'reactionformula', 'mql8[C_c] + tton[C_p]  -> 2 h[C_c] + mqn8[C_c] + so3[C_p] + tsul[C_p]','lowerBound',0,'upperBound',1000,'geneRule','');%original mql8[C_c] + tton[C_p]  -> 4 h[C_c] + mqn8[C_c] + so3[C_p] + tsul[C_p] 

model_FP = removeRxns(model_FP, {'R_SULR2pp'});
model_FP=addReaction(model_FP, 'R_SULR2pp','reactionName', 'Sulfite reductase (menaquinol) periplasmic', 'reactionformula', 'h[C_p] + 3 mql8[C_c] + so3[C_p]  -> 3 h2o[C_p] + 3 mqn8[C_c] + h2s[C_p]');%original 2 h[C_p] + 3 mql8[C_c] + so3[C_p]  -> 3 h2o[C_p] + 3 mqn8[C_c] + h2s[C_p]

model_FP = removeRxns(model_FP, {'R_SULR3_1'});
model_FP=addReaction(model_FP, 'R_SULR3_1','reactionName', 'Sulfite reductase (ferredoxin dependent)', 'reactionformula', '8 h[C_c] + so3[C_c] + 3 fdxr_42[C_c]  -> 3 h2o[C_c] + 3 fdxo_42[C_c] + h2s[C_c]');%original 8 h[C_c] + so3[C_c] + 3 fdxr_42[C_c]  -> 3 h2o[C_c] + 3 fdxo_42[C_c] + h2s[C_c]

model_FP = removeRxns(model_FP, {'R_METS'});
model_FP=addReaction(model_FP, 'R_METS','reactionName', 'Methionine synthase', 'reactionformula', '5mthf[C_c] + hcys__L[C_c] 	->	met__L[C_c] + thf[C_c]','lowerBound',0,'upperBound',1000,'geneRule','( WP_005934574_1 or WP_035394607_1 )');%original 5mthf[C_c] + hcys__L[C_c] 	->	h[C_c] + met__L[C_c] + thf[C_c]

model_FP = removeRxns(model_FP, {'R_MTHFR'});
model_FP=addReaction(model_FP, 'R_MTHFR','reactionName', '5,10-methylenetetrahydrofolate reductase (FADH2)', 'reactionformula', 'fadh2[C_c] + mlthf[C_c]  -> fad[C_c] + 5mthf[C_c]','lowerBound',0,'upperBound',1000,'geneRule','WP_005936540_1');%original fadh2[C_c] + mlthf[C_c]  -> fad[C_c] + 5mthf[C_c]

model_FP = removeRxns(model_FP, {'R_SULR3_1'});
model_FP=addReaction(model_FP, 'R_SULR3_1','reactionName', 'Sulfite reductase (ferredoxin dependent)', 'reactionformula', '7 h[C_c] + so3[C_c] + 3 fdxr_42[C_c]  -> 3 h2o[C_c] + 3 fdxo_42[C_c] + h2s[C_c]','lowerBound',0,'upperBound',1000,'geneRule','');%original 8 h[C_c] + so3[C_c] + 3 fdxr_42[C_c]  -> 3 h2o[C_c] + 3 fdxo_42[C_c] + h2s[C_c] 

model_FP = removeRxns(model_FP, {'R_APSR'});
model_FP=addReaction(model_FP, 'R_APSR','reactionName', 'APSR', 'reactionformula', 'aps[C_c] + trdrd[C_c]  -> h[C_c] + amp[C_c] + so3[C_c] + trdox[C_c]','lowerBound',0,'upperBound',1000,'geneRule','');%original aps[C_c] + trdrd[C_c]  -> 2 h[C_c] + amp[C_c] + so3[C_c] + trdox[C_c] 

model_FP = removeRxns(model_FP, {'R_ALR3'});
model_FP=addReaction(model_FP, 'R_ALR3','reactionName', 'Aldose reductase (methylglyoxal)', 'reactionformula', 'h[C_c] + nadh[C_c] + mthgxl[C_c]  -> nad[C_c] + acetol[C_c] ','lowerBound',0,'upperBound',1000,'geneRule','');% geneRule WP_005928399_1 original h[C_c] + nadph[C_c] + mthgxl[C_c]  -> nadp[C_c] + acetol[C_c] 

model_FP=addReaction(model_FP, 'R_NAt','reactionName', 'Sodium transport', 'reactionformula', 'h[C_e] + na1[C_e]  <=> h[C_c] + na1[C_c] ','lowerBound',-1000,'upperBound',1000,'geneRule','');%  

model_FP = removeRxns(model_FP, {'R_PGM_1','R_NTP1','R_LDAPAT','R_HACD3','R_KARA1i','R_HSTPTr','R_PRAMPC_1','R_IGPDH_1','R_PRAIS','R_IG3PS_1'});

model_FP = removeRxns(model_FP, {'R_NH4t','R_H2Ot','R_LLEUDr','R_ACCOAL'});%'R_ACCOAL' no GPR

% Remove duplicated reactions and loop reactions
model_FP = removeMetabolites(model_FP, '3hbycoa[C_c]');
model_FP = removeRxns(model_FP, {'R_ECOAH1_2','R_HACD1_2'});
model_FP = removeRxns(model_FP, {'R_ADK2','R_ALATA_L','R_CPK1','R_PPK2r','R_GK2'});
model_FP=addReaction(model_FP, 'R_GK2','reactionName', 'Guanylate kinase  GMPdATP', 'reactionformula', 'datp[C_c] + dgmp[C_c]  <=> dgdp[C_c] + dadp[C_c]','lowerBound',-1000,'upperBound',1000); %,'geneRule','NP_810922_1' https://enzyme.expasy.org/EC/2.7.4.8
model_FP = changeRxnBounds(model_FP,{'R_CYTK2';'R_DADK';'R_PTA2'},[0,0,-1000],'l'); %https://metacyc.org/META/NEW-IMAGE?object=EC-2.7.4.11&&redirect=T

% Reactions modified because thermodynamically infeasible
    
    model_FP = changeRxnBounds(model_FP,{'R_NDPK1';'R_NDPK2';'R_NDPK3';'R_NDPK4';'R_NDPK5';'R_NDPK6';'R_NDPK7';'R_NDPK8'},repelem(0,8),'l'); %Nucleoside-diphosphatekinase
    model_FP = changeRxnBounds(model_FP,{'R_HSDx';'R_HSDy'},repelem(0,2),'u');%L-homoserine dehydrogenase https://biocyc.org/gene?orgid=META&id=ASPKINIIHOMOSERDEHYDROGII-MONOMER
 
    
%% check for duplicated reactions (method 'FR')
checkDuplicateRxn(model_FP,'FR',1,1);% no duplicates found

%% Reactions added to the model

        % Add metabolites 
        Met={'ac[C_e]','asn__L[C_e]','asp__L[C_e]','citr__L[C_e]','cys__L[C_e]','etoh[C_e]','for[C_e]','isobuta[C_c]','isobuta[C_e]','3mb[C_c]','3mb[C_e]','ivcoa[C_c]','lac__D[C_e]','mal__L[C_e]','phe__L[C_e]','trp__L[C_e]','tyr__L[C_e]'};%'mqn8[C_e]','pydx[C_e]','ncam[C_e]'
        
        model_FP= addMultipleMetabolites(model_FP,Met);
        
        % Modification the name of Metabolites to avoid error messages later
        
        % Add Exchange reactions
        rxnToAdd.metList={'ac[C_e]','asn__L[C_e]','asp__L[C_e]','citr__L[C_e]','cys__L[C_e]','etoh[C_e]','for[C_e]','isobuta[C_e]','3mb[C_e]','lac__D[C_e]','mal__L[C_e]','phe__L[C_e]','trp__L[C_e]','tyr__L[C_e]','na1[C_e]'}; %,'mqn8[C_e]','pydx[C_e]','na1[C_e]','ncam[C_e]'
        rxnToAdd.lb    = repelem(-1000,length(rxnToAdd.metList));
        rxnToAdd.ub    = repelem(1000,length(rxnToAdd.metList));
        
        model_FP = addExchangeRxn_v2(model_FP,rxnToAdd.metList,rxnToAdd.lb,rxnToAdd.ub); %the function has been modified to add 'R_' in front of exchange reaction names

        % Add transport reactions 
            
    model_FP=addReaction(model_FP, 'R_FORt','reactionName', 'Formate diffusion', 'reactionformula', 'for[C_e] <=> for[C_c]'); 
    model_FP=addReaction(model_FP, 'R_ACtr','reactionName', 'Acetate diffusion', 'reactionformula', 'ac[C_e] <=> ac[C_c]');  
    model_FP=addReaction(model_FP, 'R_LACDt','reactionName', 'Lactate reversible transport via proton symport', 'reactionformula', 'h[C_e] + lac__D[C_e] <=> h[C_c] + lac__D[C_c]'); 
    model_FP=addReaction(model_FP, 'R_ASNt2r','reactionName', 'Asparagine reversible transport via proton symport', 'reactionformula', 'h[C_e] + asn__L[C_e] <=> h[C_c] + asn__L[C_c]');                         
    model_FP=addReaction(model_FP, 'R_ASPt2r','reactionName', 'Aspartate reversible transport via proton symport', 'reactionformula', '2 h[C_e] + asp__L[C_e] <=> 2 h[C_c] + asp__L[C_c]');             
    model_FP=addReaction(model_FP, 'R_CITRt2r','reactionName', 'L-Citrulline reversible transport via proton symport', 'reactionformula', 'h[C_e] + citr__L[C_e] <=> h[C_c] + citr__L[C_c]'); 
    model_FP=addReaction(model_FP, 'R_Cystr','reactionName', 'L-Cysteine reversible transport via proton symport', 'reactionformula', 'h[C_e] + cys__L[C_e] <=> h[C_c] + cys__L[C_c]'); 
    model_FP=addReaction(model_FP, 'R_PHEt2r','reactionName', 'L-Phenylalanine reversible transport via proton symport', 'reactionformula', 'h[C_e] + phe__L[C_e] <=> h[C_c] + phe__L[C_c]'); 
    model_FP=addReaction(model_FP, 'R_TRPt2r','reactionName', 'L-Tryptophan reversible transport via proton symport', 'reactionformula', 'h[C_e] + trp__L[C_e] <=> h[C_c] + trp__L[C_c]');
    model_FP=addReaction(model_FP, 'R_Tyrt2r','reactionName', 'L-Tyrosine reversible transport via proton symport', 'reactionformula', 'h[C_e] + tyr__L[C_e] <=> h[C_c] + tyr__L[C_c]');
    model_FP=addReaction(model_FP, 'R_Isobut','reactionName', 'Isobutyrate transport', 'reactionformula', 'isobuta[C_e] + h[C_e] <=> isobuta[C_c] + h[C_c]'); %Isobutyric acid and Isovaleric acid transporter => based on butyrate and valerate transporter
    model_FP=addReaction(model_FP, 'R_MALt','reactionName', 'Malate transport', 'reactionformula', 'mal__L[C_e] + h[C_e] -> mal__L[C_c] + h[C_c]'); %FAEPRAA2165_RS02675 from https://biocyc.org/login.html?redirect=https%3A%2F%2Fbiocyc.org%2Fgroup%3Fid%3D%3AALL-PROTEINS-6%26org-id%3DFAECPRAU model
    model_FP=addReaction(model_FP, 'R_Isoval','reactionName', 'Isovalerate transport', 'reactionformula', '3mb[C_e] + h[C_e] <=> 3mb[C_c] + h[C_c]'); %
    model_FP=addReaction(model_FP, 'R_Succ_Asp_symport','reactionName', 'Succinate transport', 'reactionformula', 'succ[C_c] + asp__L[C_e] <=> asp__L[C_c] + succ[C_e]','lowerBound',-1000,'upperBound',1000, 'geneRule', 'WP_005934526_1'); %From Uniprot and Expasy (4.1.2.4)
    model_FP=addReaction(model_FP, 'R_SUCD1','reactionName', 'Succinate synthesis', 'reactionformula', 'fum[C_c] + fadh2[C_c] -> succ[C_c] + fad[C_c]'); %From Thiele model
    model_FP=addReaction(model_FP, 'R_BTCOAACCOAT','reactionName', 'Butyryl-CoA:acetate CoA-transferase', 'reactionformula', 'btcoa[C_c] + ac[C_c] -> accoa[C_c] + but[C_c]'); %No reaction producing butyrate was described in the model
    model_FP=addReaction(model_FP, 'R_Ibutsyn','reactionName', 'Isobutyrate synthesis', 'reactionformula', '4mop[C_c] + h2o[C_c] -> isobuta[C_c] + ac[C_c]'); %Free gibbs = -14.957 kcal/mol. Described in Firmicutes/Clostridium Sporogens. (METACYC leucine degradation II).
    model_FP=addReaction(model_FP, 'R_IvalcoAsyn','reactionName', 'IsovalerateCoA synthesis', 'reactionformula', '4mop[C_c] + coa[C_c] + nad[C_c] -> ivcoa[C_c] + co2[C_c] + nadh[C_c]'); %RHEA 25177; also in Thiele I model
    model_FP=addReaction(model_FP, 'R_Ivalsyn','reactionName', 'Isovalerate synthesis', 'reactionformula', 'ivcoa[C_c] + amp[C_c] + ppi[C_c] + h[C_c] -> 3mb[C_c] + coa[C_c] + atp[C_c]'); %RHEA 46186
 
    %% Add formulas and Names of metabolites added before
    
        Names= {'Acetate','L-Asparagine','L-Aspartate','L-Citrulline','L-Cysteine','Ethanol','Formate','Isobutyrate','Isobutyrate','3-Methylbutanoic acid','3-Methylbutanoic acid','IsovalerateCoA','D-Lactate','L-Malate','L-Phenylalanine','L-Tryptophan','L-Tyrosine'};%'Menaquinone 8','Pyridoxal','Nicotinamide'
        Formulas= {'C2H3O2','C4H8N2O3','C4H6NO4','C6H13N3O3','C3H7NO2S','C2H6O','CHO2','C4H8O2','C4H8O2','C5H10O2','C5H10O2','C26H40N7O17P3S','C3H5O3','C4H4O5','C9H11NO2','C11H12N2O2','C9H11NO3'};%'C51H72O2','C8H9NO3','C6H6N2O'
    
        for i=1:numel(Met)            
            pos = strcmp(model_FP.mets,char(Met(i)));
            k= find(pos);
            model_FP.metNames(k)= Names(i);
            model_FP.metFormulas(k) = Formulas(i);
        end
        
%% Change the ATPM name for 
    model_FP = removeRxns(model_FP, {'R_ATPM'});
    model_FP=addReaction(model_FP, 'R_ngam','reactionName', 'Non-growth_ATP_maintenance', 'reactionformula', 'atp[C_c] + h2o[C_c]  -> adp[C_c] + pi[C_c] + h[C_c]','geneRule','spontaneous','lowerBound',0,'upperBound',1000);

%% Model constrains
    
            %Amino Acids
            AA= {'R_EX_citr__L_e';'R_EX_ala__L_e';'R_EX_arg__L_e';'R_EX_asn__L_e';'R_EX_asp__L_e';'R_EX_cys__L_e';'R_EX_gln__L_e';'R_EX_glu__L_e';'R_EX_gly_e';'R_EX_his__L_e';'R_EX_ile__L_e';'R_EX_leu__L_e';'R_EX_lys__L_e';'R_EX_met__L_e';'R_EX_phe__L_e';'R_EX_pro__L_e';'R_EX_ser__L_e';'R_EX_thr__L_e';'R_EX_trp__L_e';'R_EX_tyr__L_e';'R_EX_val__L_e'};
            model_FP = changeRxnBounds(model_FP,AA,[repelem(0,length(AA))],'u'); 
            model_FP = changeRxnBounds(model_FP,AA,[repelem(-1000,length(AA))],'l'); 
    
            AA_other={'R_EX_ala__D_e';'R_EX_ser__D_e'};
            model_FP = changeRxnBounds(model_FP,AA_other,[repelem(0,length(AA_other))],'u'); 
            model_FP = changeRxnBounds(model_FP,AA_other,[repelem(0,length(AA_other))],'l'); 
           
            % Nucleic Acids
            NA= {'R_EX_ade_e';'R_EX_gua_e';'R_EX_ura_e'};
            model_FP = changeRxnBounds(model_FP,NA,[repelem(0,length(NA))],'u'); 
            model_FP = changeRxnBounds(model_FP,NA,[repelem(-1000,length(NA))],'l');
            
            % Carbohydrates
            CH= {'R_EX_arab__D_e';'R_EX_arab__L_e';'R_EX_arbt_e';'R_EX_cellb_e';'R_EX_drib_e';'R_EX_fru_e';'R_EX_galur_e';'R_EX_glc__D_e';'R_EX_lcts_e';'R_EX_madg_e';'R_EX_mbdg_e';'R_EX_melib_e';'R_EX_pullulan1200_e';'R_EX_raffin_e';'R_EX_rib__D_e';'R_EX_sbt__D_e'};
            model_FP = changeRxnBounds(model_FP,CH,[repelem(0,length(CH))],'b');

            % Mucin Component (carbohydrate part)
            MCP= {'R_EX_acgam_e';'R_EX_gal_e';'R_EX_man_e';'R_EX_anhgm_e';'R_EX_acnam_e'}; %According to Thiele mannose is not in mucin 
            model_FP = changeRxnBounds(model_FP,MCP,[0,0,0,0,0],'l');%repelem(-1000,3),0,0
            model_FP = changeRxnBounds(model_FP,MCP,[repelem(0,5)],'u');

            %Ions and Vitamins
            IV={'R_EX_hco3_e';'R_EX_na1_e';'R_EX_no2_e';'R_EX_4abz_e';'R_EX_nmn_e';'R_EX_pydam_e';'R_EX_ca2_e';'R_EX_cl_e';'R_EX_cobalt2_e';'R_EX_cu2_e';'R_EX_fe2_e';'R_EX_fe3_e';'R_EX_h_e';'R_EX_hco3_e';'R_EX_k_e';'R_EX_mg2_e';'R_EX_mn2_e';'R_EX_nh4_e';'R_EX_pi_e';'R_EX_pnto__R_e';'R_EX_pydam_e';'R_EX_so4_e';'R_EX_ribflv_e';'R_EX_so4_e';'R_EX_zn2_e';'R_EX_s_e';'R_EX_mobd_e'};% cobalt2, cu2, fe2, fe3, zn2, panthotenate (pnto), pyridoxal (pydx) and riboflavin (ribfl) are not added in the medium but can be in Yeast Extract
            model_FP = changeRxnBounds(model_FP,IV,[repelem(1000,length(IV))],'u'); 
            model_FP = changeRxnBounds(model_FP,IV,[repelem(-1000,length(IV))],'l'); 
            
            model_FP = changeRxnBounds(model_FP, 'R_EX_no2_e',0,'b');%arbitrary value              
            model_FP = changeRxnBounds(model_FP, 'R_EX_s_e',0,'u');%arbitrary value  
            model_FP = changeRxnBounds(model_FP, 'R_EX_hco3_e',0,'u');%arbitrary value  

            % Gaz and water
            Gaz={'R_EX_co2_e';'R_EX_h2s_e';'R_EX_o2_e';'R_EX_h2o_e'};
            model_FP = changeRxnBounds(model_FP,Gaz,[1000,0,0,1000],'u'); 
            model_FP = changeRxnBounds(model_FP,Gaz,[0,repelem(0,2),-1000],'l'); 
          
            % Lipid and Membrane components
            LipMc={'R_EX_hxa_e';'R_EX_LalaDgluMdap_e';'R_EX_malttr_e';'R_EX_maltttr_e';'R_EX_murein4p3p_e';'R_EX_murein4p4p_e';'R_EX_murein4px4p_e';'R_EX_murein4px4px4p_e';'R_EX_murein5p3p_e';'R_EX_murein5p4p_e';'R_EX_murein5p5p_e';'R_EX_murein5p5p5p_e';'R_EX_murein5px3p_e';'R_EX_murein5px4p_e';'R_EX_murein5px4px4p_e';'R_EX_uaagmda_e';'R_EX_udcpdp_e';'R_EX_udcpo5_e';'R_EX_udcpp_e'};
            model_FP = changeRxnBounds(model_FP,LipMc,[repelem(0,length(LipMc))],'b');
                        
            % Siderophore (enterobactin)
            Sid={'R_EX_enter_e';'R_EX_feenter_e';'R_EX_salchs4_e';'R_EX_salchs4fe_e'};
            model_FP = changeRxnBounds(model_FP,Sid,[repelem(0,length(Sid))],'b');

             % Required fluxes to optimize Biomass
            model_FP = changeRxnBounds(model_FP,{'R_EX_bz_e';'R_EX_glyclt_e';'R_EX_icit_e';'R_EX_tton_e';'R_EX_coa_e'},[1,1,0,0,0],'u'); 
            model_FP = changeRxnBounds(model_FP,{'R_EX_bz_e';'R_EX_glyclt_e';'R_EX_icit_e';'R_EX_tton_e';'R_EX_coa_e'},[-1,0,0,0,0],'l'); %if bz=0 or/and tton=0 => No growth;
       
            % Unclassified
            model_FP = changeRxnBounds(model_FP,'R_EX_acald_e',0,'u'); %in the media
            model_FP = changeRxnBounds(model_FP,'R_EX_acald_e',0,'l'); 
            model_FP = changeRxnBounds(model_FP, 'R_EX_fum_e',0,'u');%arbitrary value
            model_FP = changeRxnBounds(model_FP, 'R_EX_fum_e',-1000,'l');%arbitrary value              
          
            Unclass={'R_EX_15dap_e';'R_EX_2pglyc_e';'R_EX_3oxoadp_e';'R_EX_4hphac_e';'R_EX_5mtr_e';'R_EX_acac_e';'R_EX_alaala_e';'R_EX_cgly_e';'R_EX_cyst__L_e';'R_EX_etha_e';'R_EX_hqn_e';'R_EX_hxan_e';'R_EX_id3acald_e';'R_EX_meoh_e';'R_EX_spmd_e';'R_EX_taur_e';'R_EX_urate_e';'R_EX_xan_e'};
            model_FP = changeRxnBounds(model_FP,Unclass,[repelem(0,length(Unclass))],'b');
                                
            %%  Experimental rates

            model_FP = changeRxnBounds(model_FP, 'R_EX_ac_e',1.32,'u');
            model_FP = changeRxnBounds(model_FP, 'R_EX_ac_e',0.69,'l');

            model_FP = changeRxnBounds(model_FP, 'R_EX_but_e',2.51,'u');
            model_FP = changeRxnBounds(model_FP, 'R_EX_but_e',1.68,'l');
  
            model_FP = changeRxnBounds(model_FP, 'R_EX_3mb_e',0.11,'u');  
            model_FP = changeRxnBounds(model_FP, 'R_EX_3mb_e',0,'l');  
           
            model_FP = changeRxnBounds(model_FP, 'R_EX_lac__D_e',0.06,'u');  
            model_FP = changeRxnBounds(model_FP, 'R_EX_lac__D_e',0.01,'l');  
            
            model_FP = changeRxnBounds(model_FP, 'R_EX_cit_e',0.56,'u');
            model_FP = changeRxnBounds(model_FP, 'R_EX_cit_e',-0.56,'l');

            model_FP = changeRxnBounds(model_FP, 'R_EX_succ_e',0.28,'u');
            model_FP = changeRxnBounds(model_FP, 'R_EX_succ_e',0.19,'l');
         
            model_FP = changeRxnBounds(model_FP, 'R_EX_isobuta_e',0.1,'u'); 
            model_FP = changeRxnBounds(model_FP, 'R_EX_isobuta_e',0,'l');         
          
            model_FP = changeRxnBounds(model_FP, 'R_EX_for_e',1.370,'u');
            model_FP = changeRxnBounds(model_FP, 'R_EX_for_e',0.81,'l');                  
            
            model_FP = changeRxnBounds(model_FP, 'R_EX_ppa_e',0.15,'u');
            model_FP = changeRxnBounds(model_FP, 'R_EX_ppa_e',-0.01,'l');  

            model_FP = changeRxnBounds(model_FP, 'R_EX_etoh_e',0,'u');
            model_FP = changeRxnBounds(model_FP, 'R_EX_etoh_e',0,'l');

            model_FP = changeRxnBounds(model_FP, 'R_EX_mal__L_e',-0.029,'u');
            model_FP = changeRxnBounds(model_FP, 'R_EX_mal__L_e',-0.22,'l');

         
%% Biomass reaction modification.
        % Biomass GAM modified to reach experimental growth rate
        
            model_FP = addReaction(model_FP,'Growth','reactionName', 'Biomass reaction','reactionFormula','0.000215957 10fthf[C_c] + 26 atp[C_c] + 0.1 amp[C_c] + 0.254849 glu__L[C_c] + 26 h2o[C_c] + 0.00177318 nad[C_c] + 0.000557809 coa[C_c] + 0.000432883 nadp[C_c] + 0.124 gmp[C_c] + 0.000215957 amet[C_c] + 0.436303 leu__L[C_c] + 0.233442 asp__L[C_c] + 0.000215957 fad[C_c] + 0.497466 ala__L[C_c] + 0.00504062 cl[C_c] + 0.000215957 ribflv[C_c] + 0.254849 gln__L[C_c] + 0.0886878 cys__L[C_c] + 0.208977 ser__L[C_c] + 0.245675 thr__L[C_c] + 0.0917461 his__L[C_c] + 0.000215957 thmpp[C_c] + 0.074 cmp[C_c] + 0.0136 damp[C_c] + 0.0143 dgmp[C_c] + 0.59329 gly[C_c] + 0.148832 met__L[C_c] + 0.000215957 thf[C_c] + 0.233442 asn__L[C_c] + 0.000215957 pydx5p[C_c] + 0.286451 arg__L[C_c] + 9.68419e-05 mql8[C_c] + 0.214074 pro__L[C_c] + 0.332324 lys__L[C_c] + 0.00650293 fe2[C_c] + 0.179414 phe__L[C_c] + 0.004201 so4[C_c] + 0.189029 k[C_c] + 0.281354 ile__L[C_c] + 0.000686609 cu2[C_c] + 0.00504062 ca2[C_c] + 0.00840103 mg2[C_c] + 0.000669178 mn2[C_c] + 9.68419e-05 cobalt2[C_c] + 0.000330231 zn2[C_c] + 0.00756142 fe3[C_c] + 0.0802 ump[C_c] + 0.0137 dcmp[C_c] + 0.000215957 mlthf[C_c] + 0.0135 dtmp[C_c] + 0.0550478 trp__L[C_c] + 0.409798 val__L[C_c] + 0.133541 tyr__L[C_c] + 0.0968419 uaagmda[C_c] 	->	26 adp[C_c] + 26 pi[C_c] + 26 h[C_c]')

            %% Optimization for biomass
     
    model_FP=changeObjective(model_FP,'Growth'); 
    FBAsolutionFP=optimizeCbModel(model_FP,'max','one',false); 
    printFluxVector(model_FP, FBAsolutionFP.x, 'true');
    
%% Export BT model as SMBL containing constrains
    cd '[your_path]'
    writeCbModel(model_FP,'sbml','Faecalibacterium_prausnitzii_A2_165',{'C_c','C_p','C_e'},{'cytosol','periplasm','extracellular space'}) %SBML level = 2 by default, sbmlVersion = 1 by default
    
%% Random Sampling (RAVEN)
       
        cd '[your_path]';
        model_FPrav=importModel('Faecalibacterium_prausnitzii_A2_165.xml')%120 errors

    model_FPrav = setParam(model_FPrav,'obj','Growth',1);
    sol = solveLP(model_FPrav,1);
    sol.f

    % Fix measured exchange fluxes around 5% of the value from FBA
    exIdx = getIndexes(model_FPrav,{'EX_ac_e','EX_but_e','EX_3mb_e','EX_lac__D_e','EX_cit_e','EX_succ_e','EX_isobuta_e','EX_for_e','EX_ppa_e','EX_etoh_e','EX_mal__L_e','Growth','ngam'},'rxns');%

    fluxes = sol.x(exIdx);
    model_FPrav = setParam(model_FPrav,'var',exIdx,fluxes,10);    
    [~, goodrxn] = randomSampling(model_FPrav,1,true,true,true); 
    rs = randomSampling(model_FPrav,5000,true,true,true,goodrxn); % 5000 Sampling generated
    fluxMean = full(mean(rs,2));
    fluxSD = full(std(rs,0,2));
    
    %% Save FBA predicted fluxes
    
    cd '[your_path]';
 
    out=table(model_FPrav.rxns,model_FPrav.rxnNames,constructEquations(model_FPrav),sol.x);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Flux'};
    out(out.Flux==0,:) =[];
    writetable(out,'FP_FBA.xlsx');


    %%  Save FVA predicted fluxes (mean+SD) and random sampling matrix
    clear out;out = fluxMean;
    clear out2;out2 = fluxSD;
    
    cd '[your_path]';

    out =table(model_FPrav.rxns,model_FPrav.rxnNames, constructEquations(model_FPrav), out ,out2);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Mean_RS' 'SD_RS'};
    out(out.Mean_RS==0,:) =[];
    writetable(out,'FP_RS.xlsx');

    rs(~any(rs,2),:) =[];
    RS_df=   array2table(full(rs));
    RS_df = [out(:,[1:3]),RS_df];
    writetable(RS_df,'FP_RS_5000samples.xlsx');

%% NADH and NADPH reactions
i= 'Nicotinamide adenine dinucleotide phosphate - reduced';
    [fluxes, rxnIdx] = getMetProduction(model_FPrav,i,fluxMean,true);   
    clear out;
    out.rxns    =model_FPrav.rxns(rxnIdx);
    out.rxnNames= model_FPrav.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model_FPrav,rxnIdx);
    out.fluxes  = num2cell(fluxes); 
    out = table(out.rxns,out.rxnNames,out.rxnEqns,out.fluxes);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Flux'};
    writetable(out,'FP_NADPH.xlsx');
    
i= 'Nicotinamide adenine dinucleotide - reduced';
    [fluxes, rxnIdx] = getMetProduction(model_FPrav,i,fluxMean,true);
    clear out;
    out.rxns    =model_FPrav.rxns(rxnIdx);
    out.rxnNames= model_FPrav.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model_FPrav,rxnIdx);
    out.fluxes  = num2cell(fluxes); 
    out = table(out.rxns,out.rxnNames,out.rxnEqns,out.fluxes);
    out.Properties.VariableNames = {'Reaction' 'Description' 'Equations' 'Flux'};
    writetable(out,'FP_NADH.xlsx');

%% Individual depletion analysis

    EX_RS= out(contains(out.(1),'EX_'),:);%32 rows
    temp= EX_RS(EX_RS.Mean_RS <0,:);
    IV2=strrep(IV,'R_','');
    idx = ismember(temp.Reaction,IV2); %get indexes of ions and vitamins metabolites
    temp(idx,:)=[];
    
    % Remove exp metabolites
    Exp_meta={'EX_ac_e','EX_but_e','EX_cit_e','EX_succ_e','EX_for_e','EX_etoh_e','EX_ppa_e','EX_lac__D_e','EX_isobuta_e','EX_mal__L_e','EX_3mb_e'};
    idx2 = ismember(temp.Reaction,Exp_meta); %get indexes of ions and vitamins metabolites
    temp(idx2,:)=[];%30
    
    % Remove water
    temp(strcmp(temp.Reaction,'EX_h2o_e'),:)=[]; %29
    
    temp.Reaction = strcat('R_', temp.Reaction);
    temp.name = strrep(temp.Reaction,'R_EX_','');
    
    cd '[your_path]';
    %first carbon source was threonine
    GR_prd = table();
    
    for i=1:length(temp.Reaction)        
        model_FP = changeRxnBounds(model_FP, temp.Reaction(i),0,'b');% 
        FBAsolution_temp=optimizeCbModel(model_FP,'max','one',false); % Optimize for GR
        disp([temp.Reaction(i) FBAsolution_temp.f])
        GR_prd.GR(i)= FBAsolution_temp.f;
        GR_prd.Depelted_Metabolite(i)= temp.Reaction(i);
            if isnan(FBAsolution_temp.f)~=1
                out=table(model_FP.rxns,model_FP.rxnNames,FBAsolution_temp.x);
                out.Properties.VariableNames = {'Reaction' 'Description' 'Flux'};
                out(out.Flux==0,:) =[];           
                writetable(out,string(strcat('FP_FBA_WO_',temp.Reaction(i),'.xlsx')));
            end
         model_FP = changeRxnBounds(model_FP,  temp.Reaction(i),-1000,'l');%
    end

writetable(GR_prd,string(strcat('FP_FBA_GR_indDepletion.xlsx')));
