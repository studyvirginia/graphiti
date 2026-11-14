#pragma once
#define COFIT_CBRT2 1
#define COFIT_VSP 1
#define COFIT_GJ1 1
#define COFIT_FUSE2 1
#define COFIT_FMA 1
#define COFIT_VREUSE 1
#import <Foundation/Foundation.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
namespace app_env_early {
struct Apply {
  Apply() {
    NSBundle* b = NSBundle.mainBundle;
    if (![b.bundleIdentifier isEqualToString:@"com.graphart.app"]) return;
    NSString* res = b.resourcePath;
    std::ifstream ef([res stringByAppendingPathComponent:@"in/env.txt"].UTF8String);
    std::string all, tok; std::getline(ef, all);
    std::istringstream is(all);
    while (is >> tok) { const size_t e = tok.find('='); if (e != std::string::npos) setenv(tok.substr(0, e).c_str(), tok.substr(e + 1).c_str(), 1); }
    setenv("SICIS", "1", 1); setenv("FINAL_PICK", "1", 1); setenv("SINGLE_FAST", "1", 1); setenv("S12_FREET", "1", 1);
    if (!getenv("H242D_FIX")) setenv("H242D_FIX", "0.35", 1);
    setenv("H930_MEMFIX", "1", 1);
    setenv("FJN_CISFUSE_HOST", "1", 1);
    {
      const NSInteger nm = [[NSUserDefaults standardUserDefaults] integerForKey:@"noiseModelMode"];
      if (nm == 2) { unsetenv("SICIS"); unsetenv("FJN_CISFUSE_HOST"); setenv("GRAPHITI_NOISE_OFF", "1", 1); }
    }
    if (!getenv("POST_SICIS_GPU")) setenv("POST_SICIS_GPU", "0", 1);
    if (!getenv("NF_COLL")) setenv("NF_COLL", "3", 1);
    if (!getenv("FP_KFIT")) setenv("FP_KFIT", "1", 1);
    setenv("FREET_QSPLIT", "1", 1);
    if (!getenv("H117B")) setenv("H117B", "1", 1);
    if (!getenv("SINGLE_FAST_PIN")) setenv("SINGLE_FAST_PIN", "1", 1);
    if (!getenv("FJN_CISFUSE_CHECK")) setenv("FJN_CISFUSE_CHECK", "0", 1);
    if (!getenv("FJN_COFIT")) setenv("FJN_COFIT", "300", 1);
    if (!getenv("FJN_COFIT_ETHREADS")) setenv("FJN_COFIT_ETHREADS", "4", 1);
    if (!getenv("FJN_COFIT_PQOS")) setenv("FJN_COFIT_PQOS", "2", 1);
    if (!getenv("K_QOS")) setenv("K_QOS", "3", 1);
    if (!getenv("COFIT_LAZY")) setenv("COFIT_LAZY", "1", 1);
    if (!getenv("FJN_LPT")) setenv("FJN_LPT", "1", 1);
    if (!getenv("FX_ASYNC")) setenv("FX_ASYNC", "1", 1);
    if (!getenv("SS3")) setenv("SS3", "1", 1);
    {
      static const char* const kv[][2] = {
        {"FJN_SCR_HOST", "1"}, {"FJN_FUSE_HOST", "2"}, {"FJN_WD_HOST", "1"}, {"H182_DG", "0"}, {"FJN_BINMERGE", "1"}, {"FJN_CBUCKET", "128"},
        {"NOGN", "2"}, {"NOGN_LEAN", "1"}, {"H744S", "0"}, {"FJN_HC2", "26"}, {"BU_NOCLK", "1"}, {"FJN_F8H", "1"}, {"FJN_F8C", "1"}, {"COFIT_F8H", "1"},
        {"FJN_COFIT_AUTO", "2"}, {"FJN_COFIT_PQOS", "1"}, {"FJN_COFIT_EQOS", "2"}, {"REC_PRE", "2"}, {"FJN_ARCH", "1"}, {"FJN_AGHC", "9"}, {"CAP_ASYNC", "1"}, {"RET_FAST", "1"}};
      for (const auto& e : kv) setenv(e[0], e[1], 1);
    }
    setenv("SKETCH_SIGMA_V2_BIN", [res stringByAppendingPathComponent:@"in/sigma_v2_zi5.bin"].UTF8String, 1);
    NSString* arch = [[NSSearchPathForDirectoriesInDomains(NSCachesDirectory, NSUserDomainMask, YES) firstObject] stringByAppendingPathComponent:[NSString stringWithFormat:@"psos_%llu.bin", (unsigned long long)[[[NSFileManager defaultManager] attributesOfItemAtPath:[res stringByAppendingPathComponent:@"kd/fjn_h950.metallib"] error:nil] fileSize] ^ (unsigned long long)[[[NSFileManager defaultManager] attributesOfItemAtPath:[res stringByAppendingPathComponent:@"kd/scan.metallib"] error:nil] fileSize] << 20]];
    setenv("H178F", arch.UTF8String, 1); setenv("H178W", arch.UTF8String, 1);
    setenv("H178", "1", 1);
    setenv("PSO_PAR", "1", 1);
    setenv("STROKE_IN", [NSTemporaryDirectory() stringByAppendingPathComponent:@"graphiti_stroke.txt"].UTF8String, 1);
  }
};
__attribute__((constructor(101))) static void apply_first() { static Apply a; }
}
