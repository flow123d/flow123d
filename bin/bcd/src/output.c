#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "struct.h"
#include "bcd.h"

//=============================================================================
//		PRINT FLOW BOUNDARY CONDITION FILE
//=============================================================================
void print_FBC_file(struct Problem *problem)
{
	FILE * FBC;
	struct Side *sde;
	char form1[MAXBUFF],form2[MAXBUFF];
	int i = 0;

	sprintf(form1,"%s%s%s","%d\t%d\t",problem->dbl_fmt,"\t%d\t%d\t%d\t%d\t%d\n"); //Dirichlet + Neumann
	sprintf(form2,"%s%s%s%s%s","%d\t%d\t",problem->dbl_fmt,"\t",problem->dbl_fmt,"\t%d\t%d\t%d\t%d\t%d\n"); // Newton

	FBC = open_file(problem->fbc_file,"w");
	fprintf(FBC,"$BoundaryFormat\n");
	fprintf(FBC,"1.0\t0\t8\n");
	fprintf(FBC,"$EndBoundaryFormat\n");
	fprintf(FBC,"$BoundaryConditions\n");
	fprintf(FBC,"%d\n",problem->n_wboundaries);
	FOR_BOUNDARY_SIDES(sde)
		if(sde->write == 1){
			if(sde->fbc_type != 3)
				fprintf(FBC,form1,i,sde->fbc_type,sde->fbc[0],2,sde->element->id,sde->id,1,sde->fbc_tag);
			else
				fprintf(FBC,form2,i,sde->fbc_type,sde->fbc[0],sde->fbc[1],2,sde->element->id,sde->id,1,sde->fbc_tag);
			i++;
		}
	fprintf(FBC,"$EndBoundaryConditions\n");
	fclose(FBC);
}
//=============================================================================
//		PRINT FLOW BOUNDARY CONDITION FILE
//=============================================================================
void print_TSO_file(struct Problem *problem)
{
	FILE * TSO;
	struct Element *ele;
	char form[MAXBUFF];
	int i = 0;
	int j;

	sprintf(form,"%s\t%s\t%s\t",problem->dbl_fmt,problem->dbl_fmt,problem->dbl_fmt);

	TSO = open_file(problem->tso_file,"w");
	fprintf(TSO,"$TransportSourceFormat\n");
	fprintf(TSO,"1.0\t0\t8\n");
	fprintf(TSO,"$EndTransportSourceFormat\n");
	fprintf(TSO,"$TransportSources\n");
	fprintf(TSO,"%d\n",problem->n_elm);
	FOR_ELEMENTS(ele){
				fprintf(TSO,"%d\t%d\t",ele->id,ele->id);
				for(j=0;j<problem->n_subst;j++)
					fprintf(TSO,form,0.0,0.0,0.0);
				fprintf(TSO,"\n");
			i++;
		}
	fprintf(TSO,"$EndTransportSources\n");
	fclose(TSO);
}
//=============================================================================
//		PRINT TRANSPORT BOUNDARY CONDTION FILE
//=============================================================================
void print_TBC_file(struct Problem *problem)
{
	FILE * TBC;
	struct Side *sde;
	char form[MAXBUFF];
	int i = 0;
	int j;

	sprintf(form,"%s%s","\t",problem->dbl_fmt);

	TBC = open_file(problem->tbc_file,"w");
	fprintf(TBC,"$Transport_BCDFormat\n");
	fprintf(TBC,"1.0\t0\t8\n");
	fprintf(TBC,"$EndTransport_BCDFormat\n");
	fprintf(TBC,"$Transport_BCD\n");
	fprintf(TBC,"%d\n",problem->n_wboundaries);
	FOR_BOUNDARY_SIDES(sde)
		if(sde->write == 1){
				fprintf(TBC,"%d\t%d",i,i);
				for(j=0;j<problem->n_subst;j++)
					fprintf(TBC,form,sde->tbc[j]);
				fprintf(TBC,"\n");
				i++;
		}
	fprintf(TBC,"$EndTransport_BCD\n");
	fclose(TBC);
}

//=============================================================================
//		PRINT TRANSPORT INITIAL CONDITION FILE
//=============================================================================
void print_TIC_file(struct Problem *problem)
{
	FILE * TIC;
	struct Element *ele;
	char form[MAXBUFF];
	int i = 0;
	int j;

	sprintf(form,"%s%s","\t",problem->dbl_fmt);

	TIC = open_file(problem->tic_file,"w");
	fprintf(TIC,"$ConcentrationFormat\n");
	fprintf(TIC,"1.0\t0\t8\n");
	fprintf(TIC,"$EndConcentrationFormat\n");
	fprintf(TIC,"$Concentrations\n");
	fprintf(TIC,"%d\n",problem->n_elm);
	FOR_ELEMENTS(ele){
			fprintf(TIC,"%d\t%d",ele->id,ele->id);
			for(j=0;j<problem->n_subst;j++)
				fprintf(TIC,form,ele->tic[j]);
			fprintf(TIC,"\n");
		i++;
	}
	fprintf(TIC,"$EndConcentrations\n");
	fclose(TIC);
}
//=============================================================================
//		PRINT MATERIAL FILE
//=============================================================================
void print_MTR_file(struct Problem *problem)
{
	FILE *MTR;
	struct Material *mat;
	int i = 0;


	MTR = open_file(problem->mtr_file,"w");
	fprintf(MTR,"$MaterialFormat\n");
	fprintf(MTR,"1.0\t0\t8\n");
	fprintf(MTR,"$EndMaterialFormat\n");
	fprintf(MTR,"$Materials\n");
	fprintf(MTR,"\t%d\n",problem->n_mat);
	FOR_MATERIALS(mat)
		fprintf(MTR,"\t%d\t%d%d\tx\n",mat->id,mat->dim,1);
	fprintf(MTR,"$EndMaterials\n");
	fprintf(MTR,"$Storativity\n");
	FOR_MATERIALS(mat)
			fprintf(MTR,"\t%d\t%s\n",mat->id,"s");
	fprintf(MTR,"$EndStorativity\n");
	fprintf(MTR,"$Sorption\n");
	for(i=0;i<problem->n_subst;i++)
		FOR_MATERIALS(mat)
			fprintf(MTR,"\t%d\t%d\n",mat->id,i);
	fprintf(MTR,"$EndSorption\n");
	fprintf(MTR,"$DualPorosity\n");
	FOR_MATERIALS(mat){
		fprintf(MTR,"\t%d\tn_m\tn_i",mat->id);
		for(i=0;i<problem->n_subst;i++)
			fprintf(MTR,"\ta_%d",i);
		fprintf(MTR,"\n");
	}
	fprintf(MTR,"$EndDualPorosity\n");
	fprintf(MTR,"$SorptionFraction\n");
	FOR_MATERIALS(mat){
		fprintf(MTR,"\t%d\t%s\n",mat->id,"r");
	}
	fprintf(MTR,"$EndSorptionFraction\n");
	fprintf(MTR,"$Geometry\n");
	FOR_MATERIALS(mat){
		if(mat->dim != 3)
			fprintf(MTR,"\t%d\t%d\t%s\n",mat->id,mat->dim,"1");
	}
	fprintf(MTR,"$EndGeometry\n");

	fprintf(MTR,"$Reaction\n");
	fprintf(MTR,"$EndReaction\n");

	fprintf(MTR,"$Density\t\t-UNUSED\n");
	FOR_MATERIALS(mat){
		fprintf(MTR,"\t%d\t%s\n",mat->id,"d");
	}
	fprintf(MTR,"$EndDensity\n");
	fclose(MTR);
}
