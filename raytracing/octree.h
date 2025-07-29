#pragma once
#include "triangle.h"
#include "box.h"
#include "vector.h"
#include <iostream>
#include <memory>

namespace GOAT
{
	namespace raytracing
	{
		/**
 * @brief This template class is used for internal purposes and represents an octree.
 *
 * This class is mainly used together with the surface class. When searching the intersection point of a ray with a surface described by a large number of triangles, one has
 * to test the ray against all triangles. This is often very time consuming. Therefore we used an octree algorithm. In this case the object is circumscribed by a box.
 * In a first step, the intersection between the ray and this box is tested. If the ray hits this bounding box, the box is subdivided into 8 sub-boxes (therefore the name octree).
 * The next task is to find the sub-box with the nearest intersection point. Then, this box is one again subdivided in another 8 boxes. This procedure is repeated several times
 * (according to the number MAX_RECURSIONS). Then, only the triangles on the last sub-box are tested against the ray. Depending on the total number of triangles on the surface,
 * this procedure can safe a lot of time. The Octree class is used to store the triangles in a tree-like structure, so the triangles within one box can be easily found.
 * Since this class is a template, in principal it can be used for different purposes. The Octree structure is recursive. That means each Octree is a part of the whole tree.
 * It contains up to 8 children (it can have less then 8, if some of the boxes are empty) and one parent Octree. The last Octree,  i.e. an Octree with no children,
 * is called a leaf.
 *
 */

		template <class T> class Octree
		{
		public:
			bool isLeaf; ///< Sign wether the Octree is a leaf or not. 
			Octree()
			{
				/*parent = 0;
				child = 0;*/
				nChilds = 0;
	//			Element = 0;
				nElements = 0;
				isLeaf = false;
				MAX_RECURSIONS = 3;
			}

			Octree(const Octree& O);

			~Octree()
			{
				delElements();
			}
			void delElements();///< Deletes all elements 
			void trimOctree();	 ///< Trims the tree, i.e. removing all empty childs
			void trimOctree(int rek);
			void createTree(int maxRecursions = 1, int rek = 0); ///< Prepare the tree with recursion depth maxRecursions
			void createChilds(); ///< Creates the childs of the Octree and calculates the bounding boxes
			void delChild(int i); ///< Delete the i-th child of the Octree
			void delAllChilds(); ///< Delete all children of the Octree
			/**
			 * @brief Set the bounding box of the Octree
			 * The bounding box is a circumscribing cuboid with center MP and the side lengths described by the components of d
			 * \param MP Center of the bounding box
			 * \param d Vector, where each component is the corresponding side length
			 */
			void setBoundingBox(maths::Vector<double> MP, maths::Vector<double> d);
			void setRecursiondepth(int max_recursions);
			Octree<T>* parent;

			std::vector<std::unique_ptr<Octree<T>>> child;
			std::vector<T> Element;
			size_t nElements, nChilds;
			int MAX_RECURSIONS;
			Box BBox;
		};

		template <class T> Octree<T>::Octree(const Octree& O)
		{
			if (O.parent)
				parent = O.parent;
			else
				parent = nullptr;
			nChilds = O.nChilds;
			for (const auto& ptr : O.child)
				child.push_back(std::make_unique<Octree<T>>(*ptr));

			Element = O.Element;
			nElements = Element.size();
			
			MAX_RECURSIONS = O.MAX_RECURSIONS;
			BBox = O.BBox;
			isLeaf = O.isLeaf;
		}

		template <class T> void Octree<T>::delElements()
		{
			if (nElements > 0)
				Element.clear();
			nElements = 0;
		}


		template <class T> void Octree<T>::trimOctree()
		{
			int rek = 0;
			this->trimOctree(rek);
		}

		template <class T> void Octree<T>::trimOctree(int rek)
		{
			int hnChilds = nChilds;
			bool valid[8] = { true,true,true,true,true,true, true,true };
			for (size_t i = 0; i < child.size(); i++)
			{
				
					if (child[i]->isLeaf && child[i]->Element.empty())
					{
					
						child[i].reset();
						valid[i] = false;
						nChilds--;
					}
				
				else
					child[i]->trimOctree(rek + 1);
			}

			if (hnChilds < nChilds) // Kinder wurden gel�scht ==> umsortieren notwendig
			{
				std::vector<std::unique_ptr<Octree<T> > > hChild;
				int i = 0;
				for (int j = 0; j < hnChilds; j++)
				{
					if (valid[j])
					{
						hChild.push_back(std::move(child[j]));
						i++;
					}
				}

				child = std::move(hChild);
				nChilds = child.size();
			}
		}


		template <class T> void Octree<T>::delAllChilds()
		{
			child.clear();
			nChilds = 0;
		}

		template <class T> void Octree<T>::delChild(int i)
		{
			child[i]->delAllChilds();
		}

		template <class T> void Octree<T>::createChilds()
		{
			int sx, sy, sz;
			maths::Vector<double> P, d;
			child.clear();
			for (int i = 0; i < 8; i++)
			{
				
				auto newChild = std::make_unique<Octree<T>>();
				newChild->parent = this;

				switch (i)
				{
				case 0: sx = -1; sy = -1; sz = -1; break;
				case 1: sx = +1; sy = -1; sz = -1; break;
				case 2: sx = -1; sy = +1; sz = -1; break;
				case 3: sx = +1; sy = +1; sz = -1; break;
				case 4: sx = -1; sy = -1; sz = +1; break;
				case 5: sx = +1; sy = -1; sz = +1; break;
				case 6: sx = -1; sy = +1; sz = +1; break;
				case 7: sx = +1; sy = +1; sz = +1; break;
				}
				P = BBox.P + maths::Vector<double>(sx * BBox.d[0] / 4.0, sy * BBox.d[1] / 4.0, sz * BBox.d[2] / 4.0);
				d = BBox.d / 2.0;
				newChild->BBox = Box(P, d, BBox.n);

				child.push_back(std::move(newChild));

			}
			nChilds = 8;
		}



		template <class T> void Octree<T>::createTree(int maxRecursions, int rek)
		{
			MAX_RECURSIONS = maxRecursions;
			if (rek < MAX_RECURSIONS)
			{
				createChilds();
				for (int i = 0; i < 8; i++) child[i]->createTree(maxRecursions, rek + 1);
			}
			else
			{
				nChilds = 0;
				isLeaf = true;
			}
		}

		template <class T> void Octree<T>::setBoundingBox(maths::Vector<double> MP, maths::Vector<double> d)
		{
			BBox = Box(MP, d, 1.0);
		}

		template <class T> void Octree<T>:: setRecursiondepth(int max_recursions)
		{			
			// delete the whole tree...
			delAllChilds();
			delElements();

			// ... and recreate it with the new number of
            // createTree(max_recursions);
		}



		bool checkTriangleBoxIntersection(Box B, triangle D);


		void addTriangleToTriangle(Octree<triangle>* O, triangle D, int rek);
		void addTriangleToTriangle(Octree<triangle>& O, triangle D);
		double rayOctreeIntersection(Octree<triangle>& T, maths::Vector<double> P, maths::Vector<double> k, triangle& D);
		void writeTriangleOctree(char* Fname, Octree<triangle>&);
		std::ostream& operator<< (std::ostream& os, Octree<triangle>&);
	}
}