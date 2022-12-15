#pragma once
#include "lightsrc_mc.h"
#include "lightsrc.h"
#include "objectshape.h"
#include "vector.h"
#include "detector.h"
#include "raybase.h"


namespace GOAT
{
	namespace raytracing
	{
		#define RAYTRACER_TYPE_NONE 0
		#define RAYTRACER_TYPE_PATH 1
		#define RAYTRACER_TYPE_OT   2	
		#define RAYTRACER_TYPE_PURE  3


		#define RAYTRACE_MAX_REFLEXIONS 2
		/**
		 * @brief Class defining a scene with lightsources and objects. This is a container used to inform the Raytracer about all necessary settings.
		 * Here, all informations about light sources and objects are stored. Light sources are described by classes derived from the virtual LightSrc base class. All objects are described by classes derived from the virtual
		 * class ObjectShape.
		 */
		class Scene
		{
		public:
			Scene();  ///< Standard constructor 
			Scene(const Scene& S); ///< Copy constructor
			void addObject(ObjectShape* Obj); ///< add single object to scene
			void addObjectList(int nobj, ObjectShape** obj); ///< add a list of objects to scene, nobj: number of objects
			void removeObject(int index); ///< removes object with index "index" from object list
			void removeAllObjects(); ///< removes all objects from the scene
			void addLightSource(LightSrc* ls) { addLightSource(ls, raytype); } ///< add single lightsource to scene
			void addLightSource(LightSrc* ls, int raytype); ///< add single lightsource to scene and determine ray type
			void removeLightSrc(int index); ///< remove one light source from the scene
			void removeAllLightSources(); ///< removes all light sources from the scene
			void addLightSourceRRT(LightSrc* ls, maths::Vector<std::complex<double> >Pol1, maths::Vector<std::complex<double> >Pol2);
			void addLightSourceList(int nls, LightSrc** ls); ///< add list of lightsources, nls: number of lightsources
			void addDetector(Detector* D); ///< add single detector to scene
			void addDetectorList(int nDet, Detector** D); ///< add a list of detectors to the scene, nDet: number of detectors to add
			void removeAllDetectors(); ///< remove all detectors from the scene
			void removeDetector(int index); ///< remove detector "index" from detector list
			void cleanAllDetectors(); ///< clean all detectors, i.e. all detectors are set to zero, but the detectors remain in the scene
			void setr0(double r0); ///< set the radius of the calculation space
			void setnS(std::complex<double> nS); ///< set the refractive index of the filling material in the scene
			void setnSRRT(std::complex<double> nS); ///< set the refractive index of the filling material in the scene
			void setRaytype(int raytype); ///< set the ray type for all light sources 
			void resetLS(); ///< reset all light sources. That means the counters for the rays within of the light sources are set to the first ray
			int testLS(); ///< tests, if all lightsources are outside all objects (return value: -1, if every lightsource is outside, >=0: number of the first lightsource which is inside)
			ObjectShape** Obj;   ///< List of all objects within the scene
			LightSrc** LS; ///< List of all light sources 
			LightSrc* LSRRT; ///< Light source for reversed ray tracing (RRT) 
			Detector** Det; ///< List of detectors, which are storing the electric field inside a defined area
			int nObj = 0; ///< Number of objects in the scene
			int nLS = 0;  ///< Number of light sources
			int nDet = 0; ///< Number of detectors
			std::complex<double> nS; ///< refractive index of the surrounding medium, i.e. the medium between the objects
			std::complex<double> nSRRT; ///< refractive index of the surrounding medium (RRT), i.e. the medium between the objects
			double r0=100.0; ///< Radius of the calculation space. All rays are followed within this calculation sphere.
			int raytype=LIGHTSRC_RAYTYPE_IRAY; ///< Type of the rays created by the light source. More detailed information about the available ray types and their meaning is provided 	             
		};


		/**
		* @brief This class provides all functionalities for the base raytracing code. It follows all rays from all light sources through the scene.
		*
		* The calculation space is determined by sphere with the radius r0. The Scene class is used to describe the position and types of the
		* light sources (derived from class LightSrc) and objects (derived from ObjectShape). In the moment, objects are extended volume object only. Area objects,
		* like mirrors or detectors are not implemented yet. Besides the starting and end points of each step, the electric field strength of the ray
		* at the beginning and the end of the step is stored for further calculation in derived classes. Also, the direction of the incident, the reflected and the transmitted ray is given.
		*/
		class Raytrace
		{
		public:
			Raytrace();
			Raytrace(const Scene& S);
			void setScene(const Scene& S); ///< sets Scene
			void setNumReflex(int numReflex); ///< sets number of reflexions 
			void trace(); ///< this is the starting point for the raytracing procedure
			virtual void traceEnterObject() = 0; ///< this function is called when the ray enters an object
			virtual void traceLeaveObject() = 0; ///< this function is called when the ray leaves an object
			virtual void traceStopObject() {}; ///< this function is called if no object was hidden
			void init(); ///< Initializes all values for the beginning of the Raytracing process
			int lost; ///< Rays unintentionally get lost, e.g. due to total internal reflection 
			int currentObj; ///< Number of the last object hit  (no object hit: -1)
			int currentLS; ///< Number of the current light source, which is currently in the calculation process
			maths::Vector<double> PStart, PStop; ///< Start and end point of the last step
			maths::Vector<std::complex<double> > EStart, EStop; ///< Start and end value of the electric field 
			maths::Vector<std::complex<double> > EStart2, EStop2;  ///< Start and end value of the electric field (second ray in IRay)

			///@{
			///  last hit with the objects surface: normalized direction of the incident, the reflected and the transmitted ray */
			maths::Vector<double> kin; ///< direction of the incident ray
			maths::Vector<double> kref; ///< direction of the reflected ray
			maths::Vector<double> ktrans; ///< direction of the transmitted ray
			/// @} 

			///@{ Powers stored when ray type is PRay 
			double PowRef; ///< Power of the reflected ray
			double PowIn; ///< Power of the incident ray
			double PowTrans; ///< Power of the transmitted ray
			Scene S; ///< Description of the scene
			bool useRRTParms; ///< Flag which tells the raytracing procedure if the RRT parameters of scene or the normal parameters are used within the calculation
			int type=RAYTRACER_TYPE_NONE; ///< Flag which shows which type of raytracer is selected

		private:
			/** @param ray: ray which should be traced, @param Reflexions: counter for the number of reflexions made within the ray tracing process.
				This parameter is needed to stop calculation after the maximal number of reflexions  @param recur: counter which will be set to the current recursion depth*/
			void traceOneRay(RayBase* ray, int& Reflexions, int& recur); ///< traces one ray 
			void copyRay(RayBase*& dest, RayBase* src);
			RayBase* ray; ///< current ray 
			RayBase* tray; ///< transmitted ray
			bool Abbruch; ///< flag to stop calculation
			int numReflex = RAYTRACE_MAX_REFLEXIONS;	///< current number of reflections 			
		};

		/**
		* @brief This class provides functionality to calculate the forces for optical tweezers. It is derived by the class Raytrace
		*/
		class Raytrace_OT : public Raytrace
		{
		public:
			Raytrace_OT();
			Raytrace_OT(Scene S);
			void setScene(const Scene& S);
			void trace();
			void traceLeaveObject(); ///< force calculation, when the ray leaves an object
			void traceEnterObject(); ///< force calculation, when the ray enters an object
			maths::Vector<double>* F; ///< list of the forces acting on the objects
			maths::Vector<double>* L; ///< angular momenta acting on the objects
			maths::Vector<double>** f; ///< list of the forces acting on the objects, separated for the different light sources
			maths::Vector<double>** l;///< list of the angular momenta acting on the objects, separated for the different light sources
		};


		/**
		* @brief Class which stores the start and end points of each step into a file
		*
		* This class stores all rays, which hit an object. Rays which go out of an object without hitting another one or rays which never touch an object
		* are not stored.
		*/
		class Raytrace_Path : public Raytrace
		{
		public:
			Raytrace_Path();
			Raytrace_Path(Scene); ///< Constructor with Scene Object
			void trace(std::string FName); ///< Starting point of the raytracing with the filename where the data will be stored. 
			void trace(); ///< If no Filename is given, the starting and endpoint of each ray step is stored in two arrays named P1 and P2
			void setShowOutgoingRays(bool show); ///< if true, Rays, which going out of an object without hidding a second one will be stored
			bool getShowOutgoingRays(); ///< Returns true, if Rays, which going out of an object without hidding a second one will be stored
			maths::Vector<double>* P1; ///< Here, the starting points of each step are stored if no filename is given
			maths::Vector<double>* P2; ///< Here, the ending points of each step are stored if no filename is given
			int numRays=0; ///< Number of rays stored
		private:
			void traceLeaveObject();
			void traceEnterObject();
			void storeStep();
			void traceStopObject();
			std::ofstream os;
			bool showOutgoingRays = true;
			bool storeInFile = false;

		};

		class Raytrace_pure : public Raytrace
		{
		public:
			Raytrace_pure();
			Raytrace_pure(const Scene&); ///< Constructor with Scene Object
		private:
			void traceLeaveObject() {};
			void traceEnterObject() {};
		};

	}
}