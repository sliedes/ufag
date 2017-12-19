// Ultra fast anagram generator

#include <algorithm>
#include <bitset>
#include <cassert>
#include <codecvt>
#include <fstream>
#include <iostream>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/program_options.hpp>

#include <unicode/stringoptions.h>
#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/ustream.h>
#include <unicode/utypes.h>

using std::bitset;
using std::cerr;
using std::cout;
using std::endl;
using std::forward;
using std::ifstream;
using std::make_pair;
using std::nullopt;
using std::optional;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::string;
using std::tie;
using std::unordered_map;
using std::vector;

using boost::hash_range;

namespace po = boost::program_options;

// total maximum of letters in input
constexpr int MAX_LETTERS = 128;

typedef uint8_t CharIdx;
constexpr int MAX_CHARIDX = std::numeric_limits<CharIdx>::max();

// If the given function returns false, terminate.
template <typename Fn>
bool forAllAlpha(const UnicodeString &s, Fn &&f) {
    for (int j=0, je=s.length(); j<je; j++) {
	UChar c = s[j];
	if (u_isUAlphabetic(c)) {
	    if (!f(c))
		return false;
	}
    }
    return true;
}

// We map the characters in the input to integers starting from 0 (and
// that must fit in CharIdx).
typedef unordered_map<UChar, int> CharMap;

// A multiset of characters.
class CharBag {
public:
    typedef vector<CharIdx> Vec;
    typedef bitset<MAX_CHARIDX+1> Bits;

    static optional<CharBag> fromUString(const UnicodeString &str, const CharMap &charmap);
    static optional<CharBag> fromNativeString(const string &str, const CharMap &charmap) {
	return fromUString(UnicodeString(str.c_str()), charmap);
    }
    const Vec &vec() const { return m_vec; }

    bool empty() const { return m_vec.empty(); }

    optional<CharBag> operator-(const CharBag &rhs) const;
    bool operator==(const CharBag &rhs) const;
private:
    CharBag(const Vec &vec) : m_vec(vec) {
	for (const auto &c : m_vec)
	    m_bits.set(c);
    }

    Vec m_vec;

    // A bitset of characters present in m_vec. Purely an optimization.
    Bits m_bits;
};

namespace std {
  template <> struct hash<CharBag> {
    size_t operator()(const CharBag &cs) const {
	return hash_range(cs.vec().begin(), cs.vec().end());
    }
  };
}

bool CharBag::operator==(const CharBag &rhs) const {
    return m_bits == rhs.m_bits && m_vec == rhs.m_vec;
}

optional<CharBag> CharBag::operator-(const CharBag &rhs) const {
    if ((rhs.m_bits & ~m_bits).any())
	return nullopt;
    if (rhs.m_vec.size() > m_vec.size())
	return nullopt;
    if (rhs.m_vec.size() == m_vec.size()) {
	if (rhs.m_vec == m_vec)
	    return CharBag(Vec());
	else
	    return nullopt;
    }

    Vec new_vec;
    auto lhs_it = m_vec.begin(), lhs_ie = m_vec.end();
    auto rhs_it = rhs.m_vec.begin(), rhs_ie = rhs.m_vec.end();
    while (lhs_it != lhs_ie && rhs_it != rhs_ie) {
	auto l = *lhs_it, r = *rhs_it;
	if (l < r) {
	    new_vec.emplace_back(l);
	    ++lhs_it;
	} else if (l == r) {
	    ++lhs_it;
	    ++rhs_it;
	} else
	    return nullopt;
    }

    if (lhs_it == lhs_ie && rhs_it != rhs_ie)
	return nullopt;

    while (lhs_it != lhs_ie)
	new_vec.emplace_back(*(lhs_it++));

    assert(lhs_it == lhs_ie && rhs_it == rhs_ie);
    assert(new_vec.size() == m_vec.size() - rhs.m_vec.size());
    return CharBag(new_vec);
}

optional<CharBag> CharBag::fromUString(const UnicodeString &str_, const CharMap &charmap) {
    UnicodeString str(str_);
    str.toLower();
    CharBag::Vec vec;
    if (forAllAlpha(str, [&vec, &charmap](UChar c) {
	    auto it = charmap.find(c);
	    if (it == charmap.end())
		return false;
	    vec.emplace_back(it->second);
	    return true;
	    })) {
	std::sort(vec.begin(), vec.end());
	return CharBag(vec);
    } else
	return nullopt;
}

[[maybe_unused]]
static ostream &operator<<(ostream &os, const CharBag &cs) {
    const auto &vec = cs.vec();
    os << "CharBag{";
    bool first = true;
    for (auto c : vec) {
	if (!first)
	    os << ", ";
	first = false;
	os << int(c);
    }
    return os << "}";
}

[[maybe_unused]]
static ostream &operator<<(ostream &os, const optional<CharBag> &cs) {
    if (cs)
	os << *cs;
    else
	os << "nil";
    return os;
}

static pair<vector<vector<string>>, vector<CharBag>> loadDictionary(
    const string &fname, const CharBag &cset, const CharMap &cmap) {
    ifstream f(fname);

    vector<vector<string>> words;
    vector<CharBag> charbags;
    unordered_map<CharBag, size_t> charbag_map;
    size_t count = 0;
    for (string line; std::getline(f, line);) {
	if (line.length() == 0)
	    continue;
	optional<CharBag> cs = CharBag::fromNativeString(line, cmap);
	if (cs && cset - *cs) {
	    if (cs->empty())
		continue;
	    ++count;
	    auto it = charbag_map.find(*cs);
	    if (it == charbag_map.end()) {
		words.emplace_back(vector<string>{{line}});
		charbags.emplace_back(*cs);
		charbag_map[*cs] = words.size()-1;
	    } else
		words[it->second].emplace_back(line);
	}
    }

    //cerr << "Loaded " << count << " dictionary words, " << words.size() << " distinct." << endl;
    return make_pair(words, charbags);
}

template <typename Fn>
void forAllAnagrams_iter_last(const vector<CharBag> &dict_charbags, const CharBag &charbag, Fn &&f,
			      vector<size_t> &words, size_t start_idx) {
    for (int i=start_idx, ie = dict_charbags.size(); i<ie; i++) {
	if (charbag == dict_charbags[i]) {
	    words.emplace_back(i);
	    f(words);
	    words.pop_back();
	}
    }
}

template <typename Fn>
void forAllAnagrams_iter(const vector<CharBag> &dict_charbags, const CharBag &charbag, Fn &&f, vector<size_t> &words,
			 size_t start_idx, int curr_len, int max_len) {
    if (curr_len+1 >= max_len) {
	forAllAnagrams_iter_last(dict_charbags, charbag, f, words, start_idx);
	return;
    }
    for (int i=start_idx, ie = dict_charbags.size(); i<ie; i++) {
	auto cs_opt = charbag - dict_charbags[i];
	if (cs_opt) {
	    auto cs = *cs_opt;
	    words.emplace_back(i);
	    if (cs.empty())
		f(words);
	    else
		forAllAnagrams_iter(dict_charbags, cs, forward<Fn>(f), words, i, curr_len+1, max_len);
	    words.pop_back();
	}
    }
}

template <typename Fn>
void forAllAnagrams(const vector<CharBag> &dict_charbags, const CharBag &charbag, int max_len, Fn &&f) {
    vector<size_t> words;
    forAllAnagrams_iter(dict_charbags, charbag, forward<Fn>(f), words, 0, 0, max_len);
}

// The words vector contains vectors of anagram-equivalent words.
// Output all possible combinations of them.
static void outputWords(ostream &stream, const vector<size_t> &word_idxs, const vector<vector<string>> words) {
    int size = word_idxs.size();
    vector<size_t> idxs(size);

    while (true) {
	stream << words[word_idxs[0]][idxs[0]];
	for (int i=1; i<size; i++)
	    stream << " " << words[word_idxs[i]][idxs[i]];
	stream << endl;

	// increment
	int curr = size-1;
	while (curr >= 0 && ++idxs[curr] == words[word_idxs[curr]].size()) {
	    idxs[curr] = 0;
	    curr--;
	}
	if (curr < 0)
	    break;
    }
}

static void usage(char * const *argv, po::options_description &visible) {
    cout << "Usage: " << argv[0] << " [options] sentence" << endl;
    cout << visible << endl;
    exit(0);
}

static po::variables_map parse_args(int argc, char * const *argv) {
    bool help = false;

    po::options_description visible("Allowed options"), cmdline_opt;
    visible.add_options()
	("help,h", po::bool_switch(&help)->default_value(false), "show this help")
	("dict,d", po::value<string>()->default_value("words.txt"), "dictionary (word list) to use")
	("len,l", po::value<int>()->default_value(3), "maximum anagram length in words");

    po::options_description hidden;
    hidden.add_options()
	("sentence", po::value<string>()->required(), "sentence");

    cmdline_opt.add(visible).add(hidden);

    po::positional_options_description p;
    p.add("sentence", 1);

    po::variables_map vm;
    try {
	po::store(po::command_line_parser(argc, argv).options(cmdline_opt).positional(p).run(), vm);
	po::notify(vm);
    } catch (po::error) {
	help = true;
    }

    if (help)
	usage(argv, visible);

    return vm;
}

static pair<CharMap, vector<UChar>> generateCharMap(const UnicodeString &input) {
    CharMap charmap;
    vector<UChar> reverse_charmap;

    int i = 0;
    forAllAlpha(input, [&charmap, &reverse_charmap, &i](UChar c) {
	    if (charmap.find(c) == charmap.end()) {
		if (i == MAX_CHARIDX) {
		    cerr << "Error: More than " << MAX_CHARIDX+1 << " different characters in input." << endl;
		    exit(1);
		}
		charmap[c] = i++;
		reverse_charmap.emplace_back(c);
	    }
	    return true;
	});

    return make_pair(charmap, reverse_charmap);
}

int main(int argc, char **argv) {
    std::locale::global(std::locale(""));

    auto vm = parse_args(argc, argv);

    UnicodeString input = UnicodeString(vm["sentence"].as<string>().c_str()).toLower();

    CharMap charmap;
    vector<UChar> reverse_charmap;

    tie(charmap, reverse_charmap) = generateCharMap(input);

    auto input_charbag = CharBag::fromUString(input, charmap).value();

    vector<vector<string>> dict_words;
    vector<CharBag> dict_charbags;
    tie(dict_words, dict_charbags) = loadDictionary(vm["dict"].as<string>(), input_charbag, charmap);

    forAllAnagrams(dict_charbags, input_charbag, vm["len"].as<int>(), [&dict_words](const vector<size_t> &word_idxs) {
	    assert(!word_idxs.empty());
	    outputWords(cout, word_idxs, dict_words);
	});
}
