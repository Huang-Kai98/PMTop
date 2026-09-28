// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 9/17/25.
//

#pragma once
#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;

class Config {
public:
    // 从文件加载配置
    static Config loadFromFile(const std::string& path);

    // 不带默认值（如果 key 不存在或类型不匹配，抛异常）
    template<typename T>
    T get(const std::string& key) const;

    // 获取配置值，支持 "a.b.c" 层级，带默认值
    template<typename T>
    T get(const std::string& key, const T& defaultVal) const;

private:
    explicit Config(const json& j);

    json root;
};

