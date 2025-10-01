// goat/storage/IStorageInterfaces.hpp
#pragma once

#include <string>
#include <memory>
#include "superarray.h"
#include "pulsedialog.h"

namespace GOAT {
    namespace raytracing {
        class Scene;  // Forward?Deklaration deines Datenmodells
    }

    namespace storage {

        
        /**
        * This class provides an abstract interface for data file readers
        */
        
        class IReader {
        public:
            virtual ~IReader() = default;

            /// Reads the file `filename` and returns a Scene object
            
            virtual void read(const std::string& filename) = 0;
        };

        /// Funktionszeiger-Typ für Reader?Factories (erzeugt std::unique_ptr<IReader>)
        using ReaderFactory = std::unique_ptr<IReader>(*)();

        /// Registriert eine Reader?Factory unter dem Schlüssel `name`
        void registerReader(const std::string& name, ReaderFactory factory);

        /// Erzeugt per String?Lookup den passenden Reader
        std::unique_ptr<IReader> createReader(const std::string& name);


        // ------------------------------------------------------------------------
        // Writer?Interface
        // ------------------------------------------------------------------------
        class IWriter {
        public:
            virtual ~IWriter() = default;

            /// Schreibt das Scene?Objekt in die Datei `filename`.
            virtual void write(
                const GOAT::raytracing::Scene& scene,
                GOAT::raytracing::SuperArray<GOAT::maths::Vector<std::complex<double>>>* sa,
                GOAT::visualization::VTK::pulseParameters& pp,     // gleiche cv/ref-Qualifizierung
                const std::string& filename
            ) = 0;
        };

        /// Funktionszeiger-Typ für Writer?Factories (erzeugt std::unique_ptr<IWriter>)
        using WriterFactory = std::unique_ptr<IWriter>(*)();

        /// Registriert eine Writer?Factory unter dem Schlüssel `name`
        void registerWriter(const std::string& name, WriterFactory factory);

        /// Erzeugt per String?Lookup den passende Writer
        std::unique_ptr<IWriter> createWriter(const std::string& name);

    } // namespace storage
} // namespace GOAT
