#include <flow_gtest.hh>

#include "yaml-cpp/yaml.h"

TEST(YamlCpp, parser) {
    YAML::Node config = YAML::LoadFile("test_input.yaml");
    EXPECT_TRUE(config["test_key"]);
    EXPECT_EQ(42, config["test_key"].as<int>());
}
