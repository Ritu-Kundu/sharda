#pragma once

#include <cctype>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace sharda {

inline std::string debug_query_read_text_file(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open debug artifact: " + path);
    }
    return std::string((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
}

inline size_t debug_query_find_matching(const std::string& text,
                                        size_t start,
                                        char open_char,
                                        char close_char) {
    bool in_string = false;
    bool escaped = false;
    int depth = 0;
    for (size_t index = start; index < text.size(); ++index) {
        char ch = text[index];
        if (escaped) {
            escaped = false;
            continue;
        }
        if (ch == '\\') {
            escaped = true;
            continue;
        }
        if (ch == '"') {
            in_string = !in_string;
            continue;
        }
        if (in_string) {
            continue;
        }
        if (ch == open_char) {
            depth++;
        } else if (ch == close_char) {
            depth--;
            if (depth == 0) {
                return index;
            }
        }
    }
    return std::string::npos;
}

inline std::vector<std::string> debug_query_extract_array_objects(
    const std::string& text,
    const std::string& array_key) {
    std::vector<std::string> objects;
    std::string key = "\"" + array_key + "\"";
    size_t key_pos = text.find(key);
    if (key_pos == std::string::npos) {
        return objects;
    }

    size_t array_start = text.find('[', key_pos);
    if (array_start == std::string::npos) {
        return objects;
    }
    size_t array_end = debug_query_find_matching(text, array_start, '[', ']');
    if (array_end == std::string::npos) {
        return objects;
    }

    size_t index = array_start + 1;
    while (index < array_end) {
        while (index < array_end && std::isspace(static_cast<unsigned char>(text[index]))) {
            index++;
        }
        if (index >= array_end) {
            break;
        }
        if (text[index] == ',') {
            index++;
            continue;
        }
        if (text[index] != '{') {
            index++;
            continue;
        }
        size_t object_end = debug_query_find_matching(text, index, '{', '}');
        if (object_end == std::string::npos) {
            break;
        }
        objects.push_back(text.substr(index, object_end - index + 1));
        index = object_end + 1;
    }

    return objects;
}

inline std::string debug_query_json_escape(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char ch : value) {
        switch (ch) {
        case '\\': escaped += "\\\\"; break;
        case '"': escaped += "\\\""; break;
        case '\n': escaped += "\\n"; break;
        case '\r': escaped += "\\r"; break;
        case '\t': escaped += "\\t"; break;
        default: escaped += ch; break;
        }
    }
    return escaped;
}

inline std::optional<std::string> debug_query_find_read_trace_object(
    const std::string& text,
    const std::string& read_name) {
    const std::string needle = "\"read_name\": \"" + debug_query_json_escape(read_name) + "\"";
    for (const auto& object : debug_query_extract_array_objects(text, "reads")) {
        if (object.find(needle) != std::string::npos) {
            return object;
        }
    }
    return std::nullopt;
}

inline std::optional<std::string> debug_query_find_locus_trace_object(
    const std::string& text,
    int start,
    int length) {
    const std::string start_needle = "\"local_start\": " + std::to_string(start);
    const std::string length_needle = "\"length\": " + std::to_string(length);
    for (const auto& object : debug_query_extract_array_objects(text, "loci")) {
        if (object.find(start_needle) != std::string::npos
            && object.find(length_needle) != std::string::npos) {
            return object;
        }
    }
    return std::nullopt;
}

} // namespace sharda