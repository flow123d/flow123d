#include <flow_gtest.hh>

#include "yaml-cpp/yaml.h"

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
}

