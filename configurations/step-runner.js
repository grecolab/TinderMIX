case 'TINDERMIX_LOAD_PHENOTYPE':
	step.inputs.filename.value = './test/sample/WY-14643_pheno.txt';
break;
case 'TINDERMIX_LOAD_EXPRESSION_DATA':
	step.inputs.filename.value = './test/sample/WY-14643_exp.txt';
break;
case 'TINDERMIX_COMPUTE_FC':
	step.inputs.exp_data.value = filename('TINDERMIX_LOAD_EXPRESSION_DATA', 'exp_data');
	step.inputs.pheno_data.value = filename('TINDERMIX_LOAD_PHENOTYPE', 'pheno_data');
break;
case 'TINDERMIX_CREATE_CONTOURS':
	step.inputs.fc_data.value = filename('TINDERMIX_COMPUTE_FC', 'fc_data');
	step.inputs.pheno_data.value = filename('TINDERMIX_LOAD_PHENOTYPE', 'pheno_data');
break;
case 'TINDERMIX_RUN_BMD':
	step.inputs.contour_res.value = filename('TINDERMIX_CREATE_CONTOURS', 'contour_res');
break;
case 'TINDERMIX_PLOT_N_GENE_BY_LABEL':
	step.inputs.DDRGene.value = filename('TINDERMIX_RUN_BMD', 'DDRGene');
break;
case 'TINDERMIX_PLOT_CAKE_DIAGRAM':
	step.inputs.DDRGene.value = filename('TINDERMIX_RUN_BMD', 'DDRGene');
break;
case 'TINDERMIX_CREATE_GENE_TABLE':
	step.inputs.contour_res.value = filename('TINDERMIX_CREATE_CONTOURS', 'contour_res');
	step.inputs.DDRGene.value = filename('TINDERMIX_RUN_BMD', 'DDRGene');
break;