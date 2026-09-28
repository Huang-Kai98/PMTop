// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 9/17/25.
//

#include "Config.h"
#include <fstream>
#include <stdexcept>

Config Config::loadFromFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open config file: " + path);
    }

    json j;
    in >> j;
    return Config(j);
}

Config::Config(const json& j) : root(j) {}

// 不带默认值版本（找不到 key 或类型不对 → 抛异常）
template<typename T>
T Config::get(const std::string& key) const {
    const json* cur = &root;
    size_t start = 0;
    while (start < key.size()) {
        size_t dot = key.find('.', start);
        std::string k = key.substr(start, dot - start);

        if (!cur->contains(k)) {
            throw std::runtime_error("Config key not found: " + key);
        }
        cur = &((*cur)[k]);

        if (dot == std::string::npos) break;
        start = dot + 1;
    }

    return cur->get<T>();
}

// 模板函数必须在头文件实现，或者显式实例化
template<typename T>
T Config::get(const std::string& key, const T& defaultVal) const {
    const json* cur = &root;
    size_t start = 0;
    while (start < key.size()) {
        size_t dot = key.find('.', start);
        std::string k = key.substr(start, dot - start);

        if (!cur->contains(k)) {
            return defaultVal;
        }
        cur = &((*cur)[k]);

        if (dot == std::string::npos) break;
        start = dot + 1;
    }

    try {
        return cur->get<T>();
    } catch (...) {
        return defaultVal;
    }
}

// 显式实例化常见类型
template int Config::get<int>(const std::string&, const int&) const;
template double Config::get<double>(const std::string&, const double&) const;
template bool Config::get<bool>(const std::string&, const bool&) const;
template std::string Config::get<std::string>(const std::string&, const std::string&) const;

template int Config::get<int>(const std::string&) const;
template double Config::get<double>(const std::string&) const;
template bool Config::get<bool>(const std::string&) const;
template std::string Config::get<std::string>(const std::string&) const;