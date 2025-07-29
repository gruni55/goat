#include "storageBackend.h"
#include <unordered_map>
#include <mutex>
#include <stdexcept>

namespace GOAT {
    namespace storage {

        // --------------------------------------------------
        // Meyers‑Singleton‑Accessoren
        // --------------------------------------------------
        static std::unordered_map<std::string, ReaderFactory>& readerRegistry() {
            static std::unordered_map<std::string, ReaderFactory> reg;
            return reg;
        }

        static std::unordered_map<std::string, WriterFactory>& writerRegistry() {
            static std::unordered_map<std::string, WriterFactory> reg;
            return reg;
        }

        static std::mutex& registryMutex() {
            static std::mutex m;
            return m;
        }

        // --------------------------------------------------
        // registerReader / createReader
        // --------------------------------------------------
        void registerReader(const std::string& name, ReaderFactory factory) {
            std::lock_guard<std::mutex> lock(registryMutex());
            auto& reg = readerRegistry();
            auto res = reg.emplace(name, factory);
            if (!res.second) {
                throw std::runtime_error("Reader already registered: " + name);
            }
        }

        std::unique_ptr<IReader> createReader(const std::string& name) {
            std::lock_guard<std::mutex> lock(registryMutex());
            auto& reg = readerRegistry();
            auto it = reg.find(name);
            if (it == reg.end()) {
                throw std::runtime_error("Unknown reader: " + name);
            }
            return (it->second)();
        }

        // --------------------------------------------------
        // registerWriter / createWriter
        // --------------------------------------------------
        void registerWriter(const std::string& name, WriterFactory factory) {
            std::lock_guard<std::mutex> lock(registryMutex());
            auto& reg = writerRegistry();
            auto res = reg.emplace(name, factory);
            if (!res.second) {
                throw std::runtime_error("Writer already registered: " + name);
            }
        }

        std::unique_ptr<IWriter> createWriter(const std::string& name) {
            std::lock_guard<std::mutex> lock(registryMutex());
            auto& reg = writerRegistry();
            auto it = reg.find(name);
            if (it == reg.end()) {
                throw std::runtime_error("Unknown writer: " + name);
            }
            return (it->second)();
        }

    } // namespace storage
} // namespace GOAT
