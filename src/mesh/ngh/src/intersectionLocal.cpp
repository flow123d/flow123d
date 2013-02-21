#include <iostream>
#include <algorithm>
#include "mesh/ngh/include/intersectionLocal.h"
#include "mesh/ngh/include/system.h"
#include "mesh/ngh/include/config.h"
#include "mesh/ngh/include/problem.h" //epsilon
#include <stdio.h>
#include <vector>
#include <cmath>
//#include <math.h>

int IntersectionLocal::numberInstance = 0;

int IntersectionLocal::generateId() {
    return IntersectionLocal::numberInstance++;
}

bool eps_double_equal(double a, double b) {
	return fabs(a-b) < epsilon;
}

bool IntersectionPoint::operator ==(const IntersectionPoint &IP) {
	/*cout<< "IP.coord1.begin().size: "<< IP.coord1.size() << endl;
	cout<< "el1_coord().size: "<< el1_coord().size() << endl;
	cout<< "IP.coord2.begin().size: "<< IP.coord2.size() << endl;
	cout<< "el2_coord().size: "<< el2_coord().size() << endl;
*/
	if (! equal (IP.coord1.begin(), IP.coord1.end(), el1_coord().begin(), eps_double_equal)) return false;
	if (! equal (IP.coord2.begin(), IP.coord2.end(), el2_coord().begin(), eps_double_equal)) return false;
	return true;
}

IntersectionPoint *interpolate(IntersectionPoint &A1, IntersectionPoint &A2, double t) {
    if (! (A1.el1_coord().size() == 1 && A2.el1_coord().size() == 1 ) ) {
        mythrow((char*) "Interpolation of IntersectionPoints with non line first element.", __LINE__, __FUNC__);
    }
    if (! (A1.el2_coord().size() == A2.el2_coord().size()) ) {
        mythrow((char*) "Interpolation of IntersectionPoints with non matching second element type.", __LINE__, __FUNC__);
    }
    std::vector<double> el2_coord(A1.el2_coord().size());
	for(unsigned int  i = 0; i < el2_coord.size(); ++i) {
    	el2_coord[i] = ((t - A1.el1_coord()[0]) / (A2.el1_coord()[0] - A1.el1_coord()[0])) * (A2.el2_coord()[i] - A1.el2_coord()[i]) + A1.el2_coord()[i];
	}
    return new IntersectionPoint(std::vector<double>(1,t), el2_coord);
}

IntersectionLocal::IntersectionLocal(IntersectionType i_type)
: type(i_type)
{
	id = generateId();
}

/*void IntersectionLocal::set_elements(TElement  *elem1, TElement *elem2) {
	bool pom1 = false;
	bool pom2 = false;

	// cyklus pres body pruniku => pres vsechny body v poli i_points
	for (std::vector<IntersectionPoint *>::iterator i_point = i_points.begin();
		i_point!= i_points.end();++i_point) {
			// TEST SHODY POCTU SOURADNIC S POCTEM STRAN 1.ELEMENTU
			if ((*i_point)->coord1.size() == elem1->GetNSides() - 1) { //pocet lokal. souradnic 1.elementu
				pom1 = true;
			}
			else {
				pom1 = false;
				//cout<<"Elem1 - pocet souradnic: "<< (*i_point)->coord1.size() << endl;
				//cout<<"Coord1: "<< (*i_point)->coord1[0] << endl;
				//cout<<"Coord1: "<< (*i_point)->coord1[1] << endl;
				//cout<<"Elem2 - pocet stran: "<< elem1->GetNSides() << endl;
			}
			// TEST SHODY POCTU SOURADNIC S POCTEM STRAN 2.ELEMENTU
			if ((*i_point)->coord2.size()==elem2->GetNSides() - 1) { //pocet lokal. souradnic 2.elementu
				pom2 = true;
			}
			else {
				pom2 = false;
				//cout<<"Elem1 - pocet souradnic: "<< (*i_point)->coord2.size() << endl;
				//cout<<"Coord2: "<< (*i_point)->coord2[0] << endl;
				//cout<<"Coord2: "<< (*i_point)->coord2[1] << endl;
				//cout<<"Elem2 - pocet stran: "<< elem2->GetNSides() << endl;
			}
	}
	if ((pom1==true) && (pom2==true)) {
		ele1 = elem1;
		ele2 = elem2;
	}
	else {
		mythrow((char*) "Invalid size of the vector local coordinates.", __LINE__, __FUNC__);
	}
}*/

void IntersectionLocal::add_local_coord(const std::vector<double> &coordin1, const std::vector<double> &coordin2) {
	i_points.push_back(new IntersectionPoint(coordin1, coordin2));
}

void IntersectionLocal::add_local_point(IntersectionPoint *InPoint) {
	i_points.push_back(InPoint);
	//i_points.push_back(IntersectionPoint::IntersectionPoint (InPoint));
}

void IntersectionLocal::print(FILE *out_file) {
	int size_0 = i_points.size();
	fprintf(out_file," %d", size_0); //pocet bodu pruniku N

	// cyklus pres body pruniku
	for (std::vector<IntersectionPoint *>::iterator i_point = i_points.begin();
		i_point != i_points.end();++i_point) {

			// pocet lokalnich souradnic 1. elementu a jejich vypis
	        int size_1 = (*i_point)->coord1.size();
	        fprintf(out_file, " %d", size_1);
	    	for(std::vector<double>::iterator l_coord = (*i_point)->coord1.begin();
	    		l_coord != (*i_point)->coord1.end();++l_coord) {
	    			float f_1 = *l_coord;
	    			//fprintf(out_file, " %f", f_1);
	    			fprintf(out_file, " %1.7e", f_1);
	    	}

	    	// pocet lokalnich souradnic 2. elementu a jejich vypis
	    	int size_2 = (*i_point)->coord2.size();
	    	fprintf(out_file, " %d", size_2);
		    for(std::vector<double>::iterator l_coord = (*i_point)->coord2.begin();
		    	l_coord != (*i_point)->coord2.end();++l_coord) {
		    		float f_2 = *l_coord;
		    		//fprintf(out_file, " %f", f_2);
		    		fprintf(out_file, " %1.7e", f_2);
		    }
	}
}

/*void IntersectionLocal::AddNewLocalcoord(std::vector<IntersectionPoint> &IntersectionPoint) {
	i_points.push_back(); // vlozi novy prvek na konec vektoru i_points
	i_points.pop_back();  // odeberu posledni prvek vektoru i_points
}*/

IntersectionLocal::~IntersectionLocal() {
	for (std::vector<IntersectionPoint *>::iterator i_point = i_points.begin();
		i_point != i_points.end();++i_point) {
		delete *i_point;
	}
}
