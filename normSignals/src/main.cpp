#define ARMA_64BIT_WORD
#include <armadillo>

#include <iomanip>
#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <memory>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <json/writer.h>

#include <zi/zargs/zargs.hpp>
ZiARG_string(assembly, "", "assembly");
ZiARG_string(bwFnp, "", "bigWig output name");
ZiARG_bool(debug, false, "debug");

#include "cpp/files.hpp"
#include "cpp/utility.hpp"
#include "BigWigWrapper.hpp"

namespace bib {

namespace a = arma;

struct Stats {
    double min_ = 0;
    double max_ = 0;
    double mean_ = 0;
    double stddev_ = 0;

    Stats(){}
    Stats(double mmin, double mmax, double mmean, double stddev)
        : min_(mmin), max_(mmax), mean_(mmean), stddev_(stddev)
    {}
};

class Norm{
    const bfs::path inFnp_;
    const bfs::path bwFnp_;
    const std::string assembly_;
    bfs::path chrLenFnp_;
    bfs::path blacklistFnp_;

    using BwDataT = std::map<std::string, std::vector<zentlib::Interval>>;

    uint64_t totalBases_ = 0;
    uint64_t blacklistedBases_ = 0;

    void loadData(BwDataT& bwData){
        std::cout << "loading data..." << std::endl;

        zentlib::BigWig bw(inFnp_);
        for(const auto& chrAndSize : bw.ChromsAndSizes()){
            auto chr = chrAndSize.first;
            totalBases_ += chrAndSize.second + 1;
            bwData[chr] = bw.Data(chr);
        }
    }

    auto blacklisted(){
        std::unordered_map<std::string, std::unordered_set<uint32_t>> ret;

        std::string fnp = blacklistFnp_.string();
        std::ifstream f(fnp);
        if(!f.is_open()){
            throw std::runtime_error("could not open file " + fnp + " for reading");
        }

        std::string line;
        std::string chr;
        uint32_t start;
        uint32_t end;
        while(std::getline(f, line)){
            std::istringstream ss(line);
            convert(line, ss, chr);
            convert(line, ss, start);
            convert(line, ss, end);

            for (size_t i = start; i < end; ++i) {
                ret[chr].insert(i);
            }
        }

        return ret;
    }

    void calcMeanStddev(BwDataT& bwData){
        std::cout << "loading IDR blacklist..." << std::endl;
        const auto unmappable = blacklisted();

        std::cout << "combining all values..." << std::endl;
        std::vector<float> values;
        values.reserve(totalBases_);
        for(const auto& chrAndData : bwData){
            auto& intervals = chrAndData.second;
            for(size_t i = 0; i < intervals.size(); ++i){
                const auto& d = intervals[i];
                const auto& chr = chrAndData.first;
                if(1 == unmappable.count(chr)){
                    const auto& chrBlacklist = unmappable.at(chr);
                    for (uint32_t i = d.start; i < d.end; ++i) {
                        if(unlikely(1 == chrBlacklist.count(i))){
                            ++blacklistedBases_;
                            continue;
                        }
                        values.push_back(d.val);
                    }
                } else{
                    //std::cout << "missing " << chr << std::endl;
                    for (uint32_t i = d.start; i < d.end; ++i) {
                        values.push_back(d.val);
                    }
                }
            }
        }

        writeStats(values);
    }

    Json::Value info_;
    Stats stats_;

    template <typename T>
    void showAndStore(const std::string key, const T value){
        std::cout << key << ": " << value << std::endl;
        info_[key] = value;
    }

    Stats calcAndShowAndStore(std::string prefix, const a::fvec& av){
        Stats s(a::min(av), a::max(av), a::mean(av), a::stddev(av));
        showAndStore(prefix + "-min", s.min_);
        showAndStore(prefix + "-max", s.max_);
        showAndStore(prefix + "-mean", s.mean_);
        showAndStore(prefix + "-stddev", s.stddev_);
        return s;
    }

    void writeStats(std::vector<float>& values){
        std::cout << "computing min/max/mean/stddev..." << std::endl;
        a::fvec av(values.data(), values.size(), false, true);

        calcAndShowAndStore("initial", av);

        av = a::asinh(av);
        Stats postAsinh = calcAndShowAndStore("after_asinh", av);

        av = (av - postAsinh.mean_) / postAsinh.stddev_;
        calcAndShowAndStore("after_asinh_zscore", av);

        info_["numBlacklistedBases"] = Json::UInt64{blacklistedBases_};
        info_["totalBases"] = Json::UInt64{values.size()};

        {
            bfs::path outFnp = bwFnp_.string() + ".json";
            std::ofstream out(outFnp.string());
            if(!out.is_open()){
                throw std::runtime_error("could not open file " + outFnp.string());
            }
            out << info_ << "\n";
            out.close();
            std::cout << "wrote: " << outFnp << std::endl;
        }
    }

    void normalize(BwDataT& bwData){
        // z-score normalize
        std::cout << "transforming..." << std::endl;
        for(auto& chrAndData : bwData){
            auto& intervals = chrAndData.second;
            //std::cout << "\t" << chrAndData.first << std::endl;
            #pragma omp parallel for
            for(size_t i = 0; i < intervals.size(); ++i){
                auto& d = intervals[i];
                d.val = (std::asinh(d.val) - stats_.mean_) / stats_.stddev_;
                d.val = clamp(d.val, -100.0, 100.0); // ignore extremes
            }
        }
    }

    void writeBed(const BwDataT& bwData, bfs::path bedFnp){
        std::cout << "writing transformed bed..." << std::endl;

        std::ofstream out(bedFnp.string());
        if(!out.is_open()){
            throw std::runtime_error("could not open file " +
                                     bedFnp.string());
        }
        for(auto& chrAndData : bwData){
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
        std::cout << "wrote " << bedFnp << std::endl;
    }

    void writeBigWig( bfs::path bedFnp){
        std::cout << "writing bigWig..." << std::endl;
        zentlib::BigWig::BedToBigWig(bedFnp, bwFnp_, chrLenFnp_);
        std::cout << "wrote " << bwFnp_ << std::endl;
    }

public:
    Norm(std::string inFnp, std::string bwFnp, std::string assembly)
        : inFnp_(inFnp)
        , bwFnp_(bwFnp)
        , assembly_(assembly)
    {
        bfs::path genomeFnp = "/project/umw_zhiping_weng/0_metadata/genome/";

        std::map<std::string, bfs::path> chrLenFnps{
            {"hg19", genomeFnp / "hg19.chromInfo"},
            {"mm10", genomeFnp / "mm10.chromInfo"}};
        if(0 == chrLenFnps.count(assembly)){
            throw std::runtime_error("assembly not found: '" + assembly + "'");
        }
        chrLenFnp_ = chrLenFnps.at(assembly);
        if(!bfs::exists(chrLenFnp_)){
                throw std::runtime_error("missing " + chrLenFnp_.string());
        }

        std::map<std::string, bfs::path> blacklistFnps{
                {"hg19", genomeFnp / "blacklist" / "hg19" / "wgEncodeDacMapabilityConsensusExcludable.bed"},
                {"mm10", genomeFnp / "blacklist" / "mm10" / "mm10-blacklist.bed"}};
        if(0 == blacklistFnps.count(assembly)){
                throw std::runtime_error("assembly not found: '" + assembly + "'");
        }
        blacklistFnp_ = blacklistFnps.at(assembly);
        if(!bfs::exists(blacklistFnp_)){
                throw std::runtime_error("missing " + blacklistFnp_.string());
        }

        files::ensureDir(bwFnp_);
    }

    void run(){
        if(!bfs::exists(inFnp_)){
            throw std::runtime_error("ERROR: file missing: " +
                                     inFnp_.string());
        }

        bfs::path bedFnp(bwFnp_);
        bedFnp.replace_extension(".norm.bed");

        {
            BwDataT bwData;
            loadData(bwData);
            calcMeanStddev(bwData);
            normalize(bwData);
            writeBed(bwData, bedFnp);
        }

        writeBigWig(bedFnp);
    }
};

} // namespace bib

int main(int argc, char** argv){
    zi::parse_arguments(argc, argv, true);  // modifies argc and argv
    const auto args = std::vector<std::string>(argv + 1, argv + argc);

    bib::Norm n(args.at(0), ZiARG_bwFnp, ZiARG_assembly);
    n.run();

    return 0;
}
