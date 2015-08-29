#include <flow_gtest.hh>

#include "yaml-cpp/yaml.h"

using namespace std;

TEST(YamlCpp, parser) {
    YAML::Node config = YAML::LoadFile("test_input.yaml");
    EXPECT_TRUE(config["test_key"]);
    EXPECT_EQ(42, config["test_key"].as<int>());
    EXPECT_TRUE(config["map_key"]);
    EXPECT_TRUE(config["map_key"].IsMap());
    EXPECT_TRUE(config["map_key"]["a"].IsScalar());
    EXPECT_TRUE(config["seq_key"]);
    EXPECT_TRUE(config["seq_key"].IsSequence());
    EXPECT_TRUE(config["seq_key"][1].IsScalar());
    
    EXPECT_TRUE(config["tag_key"].IsScalar());
    EXPECT_EQ("!my_int", config["tag_key"].Tag());
    EXPECT_EQ(13, config["tag_key"].as<int>());

    EXPECT_EQ(2, config["map_2"].size());
    EXPECT_EQ(0, config["map_2"]["a"].as<int>());
    EXPECT_EQ(1, config["map_2"]["b"].as<int>());
    EXPECT_EQ(3, config["map_3"].size());
    EXPECT_EQ(2, config["map_3"]["a"].as<int>());
}

TEST(YamlCpp, merge_support) {
    YAML::Node config = YAML::LoadFile("test_input.yaml");
    // YamlCpp do not have support for merge '<<:' key yet.
    // There is an pull-request from 2015, Apr 4 and there are some patches
    // with temporary solution. Until it is complete we can live without this functionality.
    EXPECT_TRUE( config["map_3"]);
    EXPECT_TRUE( config["map_3"]["b"]);
    EXPECT_TRUE( config["map_3"]["c"]);
    EXPECT_TRUE( config["map_3"]["b"].IsScalar());
    EXPECT_TRUE( config["map_3"]["c"].IsScalar());
    EXPECT_EQ("1", config["map_3"]["b"].as<string>());
    EXPECT_EQ("3", config["map_3"]["c"].as<string>());
    EXPECT_EQ(1, config["map_3"]["b"].as<int>());
    EXPECT_EQ(3, config["map_3"]["c"].as<int>());
}

