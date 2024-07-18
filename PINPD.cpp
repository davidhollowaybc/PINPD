/*
 *
 *  This file is part of the Virtual Leaf.
 *
 *  The Virtual Leaf is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The Virtual Leaf is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the Virtual Leaf.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2010 Roeland Merks.
 *
 */

#include <QObject>
#include <QtGui>

//#include "simplugin.h"

#include "parameter.h"

#include "wallbase.h"
#include "cellbase.h"
#include "PATPD.h"

#include "Pi.h"
#include "random.h"

#include "flux_function.h"


static const std::string _module_id("$Id$");

QString PATPD::ModelID(void) {
  // specify the name of your model here
  return QString( "PATPD" );
}

// return the number of chemicals your model uses
int PATPD::NChem(void) { return 3; }

// To be executed after cell division
void PATPD::OnDivide(ParentInfo *parent_info, CellBase *daughter1, CellBase *daughter2) {
  // rules to be executed after cell division go here
  // (e.g., cell differentiation rules)
}

void PATPD::SetCellColor(CellBase *c, QColor *color) { 
  // add cell coloring rules here

	double red=c->Chemical(1)/par->k[1];
	double green=c->Chemical(0)/par->k[0];
	double blue=0;
//	double blue=c->Chemical(3)/(1.+c->Chemical(3));
	color->setRgbF(red,green,blue);

}

void PATPD::CellHouseKeeping(CellBase *c) {
	// add cell behavioral rules here

	  if (c->AtBoundaryP()){

		if(c->Chemical(2) < par->rho0 && c->Chemical(2) > par->rho1){
			c->EnlargeTargetArea(par->mu);
		}
		else if(c->Chemical(2) < par->c0){
			c->EnlargeTargetArea(par->nu);
		}
		else{
		    c->EnlargeTargetArea(par->gamma);
		}

      if(c->Chemical(0) > par->vessel_inh_level && c->Chemical(2) > 0) c->EnlargeTargetArea(par->cell_expansion_rate*c->Chemical(0)/(par->vessel_expansion_rate + c->Chemical(0)));
	  
	  }
	
	if (c->Area() > 2 * c->BaseArea()) {
		c->Divide();
	} 


}

void PATPD::CelltoCellTransport(Wall *w, double *dchem_c1, double *dchem_c2) {
  // add biochemical transport rules here
  
		
	// following condition to prevent flow outside of leaf:
	if ((w->C1()->Index() >= 0) && (w->C2()->Index() >= 0) ) {
		  
	// following 2 lines to alter PAT between boundary and interior
	    double T = par->transport;
		if ((w->C2()->AtBoundaryP() && !w->C1()->AtBoundaryP()) || (w->C1()->AtBoundaryP() && !w->C2()->AtBoundaryP())) T /= par->c;
	
	/* Passive fluxes (Fick's law), rm'ed 1/25/22 

		double phi = w->Length() * par->D[0] * ( w->C2()->Chemical(0) - w->C1()->Chemical(0) );
		dchem_c1[0]+=phi;
		dchem_c2[0]-=phi; */

 // FD model of PD, 1/25/22. PD represented by Transporters_(0). Assuming Transporters1(0) = Transporters2(0); set in PDflux function.
  // cut L out of diffusion calculation, since it's now via a number of PD, 4/14/22, w->Length() *
	   
	   double PDphi = par->D[0] * w->Transporters1(0) * ( w->C2()->Chemical(0) - w->C1()->Chemical(0) );
	   dchem_c1[0]+=PDphi;
	   dchem_c2[0]-=PDphi;
	

	// PAT via PIN, from Merks 2007
	// efflux from cell 1 to cell 2
   
 	   double trans12 = ( T * w->Transporters1(1) * 
						  w->C1()->Chemical(0) / (par->ka + w->C1()->Chemical(0)) );
	
    // efflux from cell 2 to cell 1
	    double trans21 = ( T * w->Transporters2(1) * 
						  w->C2()->Chemical(0) / (par->ka + w->C2()->Chemical(0)) );
    
       dchem_c1[0] += trans21 - trans12;
       dchem_c2[0] += trans12 - trans21;

	 
	
  }
	
	// Influx at leaf "AuxinSource"
	// (as specified in initial condition). zeroth order production of auxin
	 
	if (w->AuxinSource()) { // test if wall is auxin source
       // * w->Length() rmed from following line 10/23/23
    double wall_source = par->leaf_tip_source;
	dchem_c1[0] += wall_source;
	dchem_c2[0] += wall_source;
	}  

   // Sink at leaf "AuxinSink"
	// (as specified in initial condition). Changing to 1st order decay of auxin, 8/8/17
	if (w->AuxinSink()) { // test if wall is auxin sink
    //	double wall_sink = par->sam_auxin_breakdown;
	// * w->Length();
	//if(w->C1()->Chemical(0) > wall_sink) dchem_c1[0] -= wall_sink; Old zeroth order removal, pre 8/8/17
	//if(w->C2()->Chemical(0) > wall_sink) dchem_c2[0] -= wall_sink;

	// 10/30/17: changed to biological 'sink', a cell into which auxin flows, via UTG. Therefore added zeroth order production, 
		// as well as first-order decay, to keep steady state. Not keeping 'length' as part of source term. 
		  dchem_c1[0] += par->leaf_tip_source - w->C1()->Chemical(0) * par->sam_auxin_breakdown;
		  dchem_c2[0] += par->leaf_tip_source - w->C2()->Chemical(0) * par->sam_auxin_breakdown;
		  dchem_c1[1] -= par->pin_breakdown_internal * w->C1()->Chemical(1);
		  dchem_c2[1] -= par->pin_breakdown_internal * w->C2()->Chemical(1);
	}  


}
void PATPD::WallDynamics(Wall *w, double *dw1, double *dw2) {
  // add biochemical networks for reactions occuring at walls here

	
	dw1[2] = 0.; dw2[2] = 0.; // chemical 2 unused in walls

// following block calculates flux for PD change

    double pd_net = 0.;

  	// The following block needed to calculate flux due to PIN in a time step
	
			
	double	pin_atwall_this = w->Transporters1(1);
	double	pin_atwall_adj = w->Transporters2(1);
	
	
	// The following block calculates abs_phi_tot, total flux, in order to have a current value for the wall allocation step. 
	// removed L dependence, 4/14/22, w->Length()
		
    if ((w->C1()->Index() >= 0) && (w->C2()->Index() >= 0))  {

        // Following block added back 10/31/23, combining PDphi and PATphi

        double T = par->transport;
		if ((w->C2()->AtBoundaryP() && !w->C1()->AtBoundaryP()) || (w->C1()->AtBoundaryP() && !w->C2()->AtBoundaryP())) T /= par->c;

		double trans_thisadj = (T * pin_atwall_this * w->C1()->Chemical(0) / (par->ka + w->C1()->Chemical(0)));

		double trans_adjthis = (T * pin_atwall_adj * w->C2()->Chemical(0) / (par->ka + w->C2()->Chemical(0)));

        double abs_phi_tot = abs((par->D[0] * w->Transporters1(0) * (w->C1()->Chemical(0) - w->C2()->Chemical(0))) + trans_thisadj - trans_adjthis );

        double PDphi = abs(par->D[0] * w->Transporters1(0) * ( w->C2()->Chemical(0) - w->C1()->Chemical(0) ));
	

        /* Facilitated diffusion via PD
               4/26/22, ceiling on Dij, next line */
            // D-only 'PDphi', 11/22/23
        // T + D 'abs_phi_tot', 11/23/23 */


            if(w->Transporters1(0) < par->k[2])
			{
                 pd_net = par->van3prod * pow(abs_phi_tot, 2) + par->van3autokat * pow(abs_phi_tot, 1)
							 + par->van3sat + par->sam_auxin * ((w->C1()->Chemical(0) + w->C2()->Chemical(0))/2) 
									- par->k2van3 * w->Transporters1(0);
									
				if(w->C1()->AtBoundaryP() && w->C2()->AtBoundaryP()){     
					if(w->C1()->CellType() != w->C2()->CellType()){
						pd_net = 0.;
			        }
		        }
			}

   
	dw1[0] = pd_net;
    dw2[0] = dw1[0]; 
   

   dw1[1] = PINflux(w->C1(),w->C2(),w);
   dw2[1] = PINflux(w->C2(),w->C1(),w); 
  
  }
}

void PATPD::CellDynamics(CellBase *c, double *dchem) { 
  // add biochemical networks for intracellular reactions here

// 1st order auxin production in boundary cells    	
	if(c->AtBoundaryP()){
	   if(c->Chemical(2) > par->f){
		 dchem[0] += par->aux1prod * c->Chemical(2) - par->aux1decay * c->Chemical(0);
		// dchem[1] += par->pin_prod * c->Chemical(0) - par->pin_breakdown * c->Chemical(1) - SumFluxFromWalls( c, PATPD::PINflux );
	   }

	   else{
		 dchem[0] += (- par->aux1decay * c->Chemical(0)); 
		 dchem[1] += 0.;
	   }

	   dchem[2] += par->aux_cons;
	}

	else{
		dchem[0] += (par->aux1prodmeso - par->aux1decay * c->Chemical(0));
		
	}

	dchem[1] += par->pin_prod * c->Chemical(0) - par->pin_breakdown * c->Chemical(1) - SumFluxFromWalls( c, PATPD::PINflux );
		 
}

 double PATPD::PINflux(CellBase *this_cell,
	CellBase *adjacent_cell, Wall *w) { 

   
	double wtf_add = 0.;
	double utg_add = 0.;
	double pin_net = 0.;

  	

 /*  double pin_atwall; // pick the correct side of the Wall
	   if (w->C1() == this_cell) pin_atwall = w->Transporters1(1);
	   else pin_atwall=w->Transporters2(1); */

	double pin_atwall_this, pin_atwall_adj; // pick the correct side of the Wall
	if (w->C1() == this_cell){
		pin_atwall_this = w->Transporters1(1);
		pin_atwall_adj = w->Transporters2(1);
	}
	else {
		pin_atwall_this=w->Transporters2(1);
		pin_atwall_adj=w->Transporters1(1);
	}

	 // calculate PIN translocation rate from cell to membrane; WTF (R-L and P, 05), but retaining saturated T flow of Merks07 eq1. 
	  // This calculates phi_tot, total flux, in order to have a current value for the wall allocation step. 
	// fragment of below condition, && w->C1()->CellType() == 0 && w->C2()->CellType() == 0)
	
	if ((w->C1()->Index() >= 0) && (w->C2()->Index() >= 0))  {
		
		double T = par->transport;
		if ((w->C2()->AtBoundaryP() && !w->C1()->AtBoundaryP()) || (w->C1()->AtBoundaryP() && !w->C2()->AtBoundaryP())) T /= par->c;

		double trans_thisadj = (T * pin_atwall_this * this_cell->Chemical(0) / (par->ka + this_cell->Chemical(0)));

		double trans_adjthis = (T * pin_atwall_adj * adjacent_cell->Chemical(0) / (par->ka + adjacent_cell->Chemical(0)));

        double phi_tot = (par->D[0] * w->Transporters1(0) * (this_cell->Chemical(0) - adjacent_cell->Chemical(0))) + trans_thisadj - trans_adjthis;

       // 10/24/23, took out this phiT-only
        // double phi_tot = trans_thisadj - trans_adjthis;


		if (phi_tot > 0){
//			if ((this_cell->AtBoundaryP() && this_cell->Chemical(0) > par->f) || (!this_cell->AtBoundaryP())) { USES OLD par->f, pre 11/23/17
			// next line added 2-8-18. par->f threshold for wtf add rm-ed 5/10/18, but it's been set arbitrarily high prior to this. 
				wtf_add = (this_cell->Chemical(1) / (par->kap + this_cell->Chemical(1))) * (par->e * pow(phi_tot, 2) + par->d * pow(phi_tot, 1));
		}
	}
	  	
		// version from R-L + P 05, though they had a hard stop of Pi<4 (no ij in their model)
		//	pin_flux = par->k1 * pow(phi_tot, 2) + par->d - par->k2 * pin_atwall1; 
	    
	// UTG allocation, added back in, 10/16/17. AtBoundary implements only for boundary cells. 
	
      double receptor_level = adjacent_cell->Chemical(0) * par->r / (par->kr + adjacent_cell->Chemical(0));  
	      
	     utg_add = par->k1 * this_cell->Chemical(1) * receptor_level / ( par->km + this_cell->Chemical(1) );

	// total allocation

		   if(this_cell->AtBoundaryP() && adjacent_cell->AtBoundaryP()){ 
			   if(this_cell->CellType() == adjacent_cell->CellType()){
					pin_net = utg_add + wtf_add - par->k2 * pin_atwall_this;
		       }
			   else{
				   pin_net = 0.;
			   }
		   }


           else {
               pin_net = utg_add + wtf_add + par->k[3] - par->k2 * pin_atwall_this;
		   }
	
			  
  return pin_net;  


	
}


//Q_EXPORT_PLUGIN2(PATPD, PATPD)
