#include "octree.h"
#include "triangle_box_intersection.h"
#include <fstream> 
#include <limits>
namespace GOAT
{
	namespace raytracing
	{
		bool checkTriangleBoxIntersection(Box B, triangle D)
		{
			return triangleAABBIntersection(B, D);
		}


		void addTriangleToTriangle(Octree<triangle>* O, triangle D, int rek)
		{
			maths::Vector<double> P, d;
			 
			if (rek < O->MAX_RECURSIONS) // maximale Rekursionstiefe nicht erreicht 
			{
				if (checkTriangleBoxIntersection(O->BBox, D))
				{

					O->isLeaf = false;

					for (int i = 0; i < 8; i++) addTriangleToTriangle((Octree<triangle> *) O->child[i], D, rek + 1);
				}
			}
			else // maximale Rekursionstiefe erreicht => Dreieck hinzufügen
			{
				O->isLeaf = true;
				if (checkTriangleBoxIntersection(O->BBox, D))
				{
					if (O->nElements == 0) O->Element = (triangle*)malloc(sizeof(triangle));
					else O->Element = (triangle*)realloc(O->Element, sizeof(triangle) * (O->nElements + 1));
					O->Element[O->nElements] = D;					
					O->nElements++;
				}
			}
		}

		void addTriangleToTriangle(Octree<triangle>& O, triangle D)

		{
			int rek = 0;
			addTriangleToTriangle(&O, D, rek);
		}

		double rayOctreeIntersection(Octree<triangle>& T, maths::Vector<double> P, maths::Vector<double> k, triangle& D)
		{
			double dmin = std::numeric_limits<double>::infinity();
			if (T.isLeaf || T.nChilds == 0) // Jetzt die Dreiecke testen
			{
				int hit = -1;
				double d;

				for (int i = 0; i < T.nElements; i++)
				{
					d = T.Element[i].distance(P, k);
					if (d < dmin)
					{
						hit = i;
						dmin = d;
					}
				}

				if (hit >= 0)
				{
					D = T.Element[hit];
					return dmin;
				}
				else return -1;
			}
			else
			{
				triangle dhilf;
				double d;
				maths::Vector<double> philf;
				bool found = false;
				if (T.BBox.next(P, k, philf))
				{
					for (int i = 0; i < T.nChilds; i++)
					{
						d = rayOctreeIntersection(*T.child[i], P, k, dhilf);
						if ((d > 0) && (d < dmin))
						{
							dmin = d;
							D = dhilf;
							found = true;
						}
					}
				}
				if (!found) return -1;
			}
			return dmin;
		}


		std::ostream& operator<< (std::ostream& os, Octree<triangle>& O)
		{
			if (O.isLeaf)
			{
				os << "% LEAF" << std::endl;
				os << "% " << O.nElements << std::endl;
				os << "% Bounding Box: " << O.BBox << std::endl << std::flush;
				for (int i = 0; i < O.nElements; i++)
					os << O.Element[i] << std::endl;
			}
			else
			{
				os << "% KNOT" << std::endl;
				os << "% " << O.nElements << std::endl;
				os << "% Bounding Box: " << O.BBox << std::endl << std::flush;
				if (O.nElements > 0)
					for (int i = 0; i < O.nElements; i++)
						os << O.Element[i] << std::endl << std::flush;

				if (O.nChilds > 0)
					for (int i = 0; i < O.nChilds; i++)
						os << *O.child[i] << std::endl << std::flush;
			}
			return os;
		}

		void writeOctreeItem(std::ofstream& os, Octree<triangle>& O, int Index)
		{
			os << "% Index " << Index << std::endl;
			os << "% Bounding Box: " << O.BBox << std::endl;
			os << O << std::endl;
			if (!O.isLeaf)
			{
				for (int i = 0; i < O.nChilds; i++)
					if (O.child[i] != 0)
					{
						int index = Index * 10 + i + 1;
						writeOctreeItem(os, *O.child[i], index);
					}
			}
		}

		void writeTriangleOctree(char* Fname, Octree<triangle>& O)
		{
			std::ofstream os;
			os.open(Fname);
			writeOctreeItem(os, O, 1);
			os.close();
		}
	}
}