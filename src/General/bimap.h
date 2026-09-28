//
// Created by huangkai on 25-1-20.
//

#ifndef BIMAP_H
#define BIMAP_H
#include<unordered_map>
#include <stdexcept>
#include <iostream>

template <typename KeyType, typename ValueType> class BiMap {
public:
  // 插入键值对
  void insert(const KeyType &key, const ValueType &value) {
    if (keyToValue.count(key) || valueToKey.count(value)) {
      throw std::invalid_argument("Key or value already exists in BiMap.");
    }
    keyToValue[key] = value;
    valueToKey[value] = key;
  }

  // 通过键查找值
  ValueType getValue(const KeyType &key) const {
    auto it = keyToValue.find(key);
    if (it == keyToValue.end()) {
      throw std::out_of_range("Key not found.");
    }
    return it->second;
  }

  // 通过值查找键
  KeyType getKey(const ValueType &value) const {
    auto it = valueToKey.find(value);
    if (it == valueToKey.end()) {
      throw std::out_of_range("Value not found.");
    }
    return it->second;
  }

  // 删除键值对
  void eraseByKey(const KeyType &key) {
    auto it = keyToValue.find(key);
    if (it != keyToValue.end()) {
      ValueType value = it->second;
      keyToValue.erase(it);
      valueToKey.erase(value);
    }
  }

  void eraseByValue(const ValueType &value) {
    auto it = valueToKey.find(value);
    if (it != valueToKey.end()) {
      KeyType key = it->second;
      valueToKey.erase(it);
      keyToValue.erase(key);
    }
  }

  // 更新键值对：通过键更新值
  void updateByKey(const KeyType &key, const ValueType &newValue) {
    auto it = keyToValue.find(key);
    if (it == keyToValue.end()) {
      throw std::out_of_range("Key not found.");
    }
    ValueType oldValue = it->second;
    keyToValue[key] = newValue;
    valueToKey.erase(oldValue);
    valueToKey[newValue] = key;
  }

  // 更新键值对：通过值更新键
  void updateByValue(const ValueType &value, const KeyType &newKey) {
    auto it = valueToKey.find(value);
    if (it == valueToKey.end()) {
      throw std::out_of_range("Value not found.");
    }
    KeyType oldKey = it->second;
    valueToKey[value] = newKey;
    keyToValue.erase(oldKey);
    keyToValue[newKey] = value;
  }

  // 检查键或值是否存在
  bool containsKey(const KeyType &key) const {
    return keyToValue.count(key) > 0;
  }

  bool containsValue(const ValueType &value) const {
    return valueToKey.count(value) > 0;
  }

  // 输出所有键值对
  void output() const {
    std::cout << "BiMap contents:\n";
    for (const auto &pair : keyToValue) {
      std::cout << "Key: " << pair.first << ", Value: " << pair.second << "\n";
    }
  }

private:
  std::unordered_map<KeyType, ValueType> keyToValue;
  std::unordered_map<ValueType, KeyType> valueToKey;
};



#endif //BIMAP_H
