#define ARMA_64BIT_WORD
#include <armadillo>

#include <iomanip>
#include <iostream>
#include <string>
#include <map>
#include <unordered_set>
#include <memory>
#include <unordered_map>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <json/writer.h>

#include <zi/zargs/zargs.hpp>
ZiARG_string(assembly, "", "assembly");
ZiARG_bool(debug, false, "debug");

#include "cpp/files.hpp"
#include "cpp/utility.hpp"
#include "BigWigWrapper.hpp"

namespace bib {

namespace a = arma;

class Norm{
    const bfs::path inFnp_;
    const std::string assembly_;
    bfs::path chrLenFnp_;

    std::map<std::string, std::vector<zentlib::Interval>> bwData_;
    uint64_t totalBases_ = 0;
    float mean_ = 0;
    float stddev_ = 0;

    void loadData(){
        std::cout << "loading data..." << std::endl;

        zentlib::BigWig bw(inFnp_);
        for(const auto& chrAndSize : bw.ChromsAndSizes()){
            auto chr = chrAndSize.first;
            totalBases_ += chrAndSize.second + 1;
            bwData_[chr] = bw.Data(chr);
        }
    }

    auto calcMeanStddev(){
        std::cout << "combining all values..." << std::endl;

        std::vector<float> values;
        values.reserve(totalBases_);
        for(const auto& chrAndData : bwData_){
            auto& intervals = chrAndData.second;
            for(size_t i = 0; i < intervals.size(); ++i){
                const auto& d = intervals[i];
                for (uint32_t i = d.start; i < d.end; ++i) {
                    if(0){
                        // check unmappable....
                    }
                    values.push_back(d.val);
                }
            }
            //std::cout << "\t" << chrAndData.first << " " << values.size() << "\n";
        }

        std::cout << "computing mean/stddev..." << std::endl;
        a::fvec av(values.data(), values.size(), false, true);
        mean_ = a::mean(av);
        stddev_ = a::stddev(av);
    }

    void normalize(){
        // z-score normalize
        std::cout << "transforming..." << std::endl;
        for(auto& chrAndData : bwData_){
            auto& intervals = chrAndData.second;
            //std::cout << "\t" << chrAndData.first << std::endl;
            #pragma omp parallel for
            for(size_t i = 0; i < intervals.size(); ++i){
                auto& d = intervals[i];
                d.val = (d.val - mean_) / stddev_;
                d.val = clamp(d.val, -1.96, 1.96); // ignore extremes
                d.val += 1.96; // shift up
            }
        }
    }

    void writeBed(bfs::path outFnp){
        std::cout << "writing transformed bed..." << std::endl;

        std::ofstream out(outFnp.string());
        if(!out.is_open()){
            throw std::runtime_error("could not open file " +
                                     outFnp.string());
        }
        for(auto& chrAndData : bwData_){
            //std::cout << "\t" << chrAndData.first << std::endl;
            for(auto& d : chrAndData.second){
                if(0 == d.val){
                    continue;
                }
                out << chrAndData.first << "\t"
                    << d.start << "\t"
                    << d.end << "\t"
                    << d.val << "\n";
            }
        }
        out.close();
        std::cout << "wrote " << outFnp << std::endl;
    }

    void writeBigWig(const bfs::path outFnp){
        bfs::path bwFnp = outFnp;
        bwFnp.replace_extension(".bigWig");
        std::cout << "writing bigWig..." << std::endl;
        zentlib::BigWig::BedToBigWig(outFnp, bwFnp, chrLenFnp_);
        std::cout << "wrote " << bwFnp << std::endl;
    }

public:
    Norm(std::string inFnp, std::string assembly)
        : inFnp_(inFnp)
        , assembly_(assembly)
    {
        //TODO:
        //chrLenFnp_ =
    }

    void run(){
        if(!bfs::exists(inFnp_)){
            throw std::runtime_error("ERROR: file missing: " +
                                     inFnp_.string());
        }

        loadData();
        calcMeanStddev();
        normalize();

        bfs::path outFnp = str::replace(inFnp_.string(),
                                        "encode/data/", "encode/norm");
        outFnp.replace_extension(".norm.bed");
        files::ensureDir(outFnp);

        writeBed(outFnp);
        writeBigWig(outFnp);
    }
};

} // namespace bib

int main(int argc, char** argv){
    zi::parse_arguments(argc, argv, true);  // modifies argc and argv
    const auto args = std::vector<std::string>(argv + 1, argv + argc);

    bib::Norm n(args.at(0), ZiARG_assembly);
    n.run();

    return 0;
}
