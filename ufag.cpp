// Ultra fast anagram generator

#include <algorithm>
#include <bitset>
#include <cassert>
#include <codecvt>
#include <fstream>
#include <ios>
#include <iostream>
#include <limits>
#include <locale>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gmpxx.h>

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
using std::ios;
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

typedef unsigned __int128 bigint;

namespace po = boost::program_options;

// total maximum of letters in input
constexpr int MAX_LETTERS = 128;

typedef uint8_t CharIdx;
constexpr int MAX_CHARIDX = std::numeric_limits<CharIdx>::max();

// first 168 primes
constexpr int PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
    919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997
};

// If the given function returns false, terminate.
template <typename Fn>
bool forAllAlpha(const UnicodeString &s, Fn &&f) {
    for (int j = 0, je = s.length(); j<je; j++) {
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

static size_t hash_bigint(const bigint &m) {
    // Let's decide the lowest bits are good enough, even though it's
    // guaranteed we have count(most_common_letter) trailing zero
    // bits.
    return static_cast<size_t>(m);
}

// A multiset of characters.
//
// The allowable characters are passed in the CharMap, which maps
// unicode characters to integers starting from 0.
//
// To facilitate the required multiset operations, namely subtraction
// and subset testing, we map a multiset to an integer (which may not
// fit in a fixed size integer type)
//
//    product(nth_prime(i)**char_count(i)).
//
// Now,
//
// * x is a subset of y iff x.num % y.num == 0
// * x-y maps to integer division.
//
// The order of characters is chosen such that the most common
// characters in the input map to smallest primes, which results in
// smallest numbers.
class CharBag {
public:
    static optional<CharBag> fromUString(const UnicodeString &str, const CharMap &charmap);
    static optional<CharBag> fromLowerUString(const UnicodeString &str, const CharMap &charmap);
    static optional<CharBag> fromNativeString(const string &str, const CharMap &charmap) {
	return fromUString(UnicodeString(str.c_str()), charmap);
    }
    const bigint &num() const { return m_num; }
    const int size() const { return m_size; }
    size_t hash() const { return m_hash; }

    bool empty() const { return m_size == 0; }

    bool isSubsetOf(const CharBag &other) const {
	auto cs_opt = other-*this;
	return (bool)cs_opt;
    }

    optional<CharBag> operator-(const CharBag &rhs) const;
    bool operator==(const CharBag &rhs) const;
private:
    CharBag(bigint num, int size) :
	m_num(num), m_size(size), m_hash(compute_hash(m_size, m_num)) {}

    static size_t compute_hash(const int &size, const bigint &num) {
	// As a hack to find words faster, we want longer words to map
	// to smaller hashes. Do this by mapping negated size to the
	// high bits of the hash.
	constexpr int bits = std::numeric_limits<size_t>::digits;
	size_t hi = (~static_cast<size_t>(size)) << (bits-7);
	size_t lo = hash_bigint(num) & ((size_t(1) << (bits-7)) - 1);
	return hi | lo;
    }

    bigint m_num;
    int m_size;
    size_t m_hash;
};

namespace std {
  template <> struct hash<CharBag> {
    size_t operator()(const CharBag &cs) const {
	return cs.hash();
    }
  };
}

bool CharBag::operator==(const CharBag &rhs) const {
    return m_hash == rhs.m_hash && m_num == rhs.m_num;
}

optional<CharBag> CharBag::operator-(const CharBag &rhs) const {
    if (rhs.m_size > m_size)
	return nullopt;

    if (m_size == rhs.m_size && m_num == rhs.m_num)
	return CharBag{bigint(1), 0};

    // if (!mpz_divisible_p(m_num.get_mpz_t(), rhs.m_num.get_mpz_t()))
    // 	return nullopt;
    if (m_num % rhs.m_num != 0)
	return nullopt;

    return CharBag{m_num/rhs.m_num, m_size-rhs.m_size};
}

optional<CharBag> CharBag::fromLowerUString(const UnicodeString &str, const CharMap &charmap) {
    bigint n(1);
    int size = 0;

    static_assert(sizeof(PRIMES)/sizeof(PRIMES[0]) >= MAX_LETTERS);

    if (forAllAlpha(str, [&n, &size, &charmap](UChar c) {
	    auto it = charmap.find(c);
	    if (it == charmap.end())
		return false;
	    int idx = it->second;
	    assert(idx < MAX_LETTERS);
	    n *= PRIMES[idx];
	    ++size;
	    return true;
	    })) {
	return CharBag(n, size);
    } else
	return nullopt;
}

optional<CharBag> CharBag::fromUString(const UnicodeString &str_, const CharMap &charmap) {
    UnicodeString str(str_);
    str.toLower();
    return fromLowerUString(str, charmap);
}

// [[maybe_unused]]
// static ostream &operator<<(ostream &os, const CharBag &cs) {
//     return os << "CharBag{" << cs.size() << ", " << cs.num() << "}";
// }

// [[maybe_unused]]
// static ostream &operator<<(ostream &os, const optional<CharBag> &cs) {
//     if (cs)
// 	os << *cs;
//     else
// 	os << "nil";
//     return os;
// }

// from https://stackoverflow.com/questions/17074324/
template <typename T>
void apply_permutation_in_place(std::vector<T>& vec,
				const std::vector<std::size_t>& p) {
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (done[i])
            continue;
        done[i] = true;
        std::size_t prev_j = i;
        std::size_t j = p[i];
        while (i != j) {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j = j;
            j = p[j];
        }
    }
}

static vector<char> slurp(const string &fileName) {
    ifstream ifs(fileName.c_str(), ios::binary | ios::ate);

    ifstream::pos_type size = ifs.tellg();
    ifs.seekg(0, ios::beg);

    vector<char> bytes(size);
    ifs.read(bytes.data(), size);

    return bytes;
}

// The dictionary words are sorted by length (number of alphabetic
// characters). This is important, it's used in
// forAllAnagrams_iter_last().
static pair<vector<vector<string>>, vector<CharBag>> loadDictionary(
    const string &fname, const CharBag &cset, const CharMap &cmap) {
    const auto raw_contents = slurp(fname);
    UnicodeString contents = UnicodeString(raw_contents.data(), raw_contents.size());
    contents.toLower();
    size_t contents_len = contents.length();

    vector<vector<string>> words;
    vector<CharBag> charbags;
    unordered_map<CharBag, size_t> charbag_map;

    size_t count = 0;

    size_t line_start = 0;
    while (true) {
	auto newline = contents.indexOf(UChar('\n'), line_start);

	size_t endpos;
	if (newline == -1)
	    endpos = contents_len;
	else
	    endpos = newline;

	auto str_len = endpos - line_start;
	auto line = contents.tempSubString(line_start, str_len);
	if (str_len > 0) {
	    optional<CharBag> cs = CharBag::fromLowerUString(line, cmap);
	    if (cs && cset - *cs) {
		if (cs->empty())
		    continue;
		++count;

		ostringstream line_sstream;
		line_sstream << line;

		auto it = charbag_map.find(*cs);
		if (it == charbag_map.end()) {
		    words.emplace_back(vector<string>{{line_sstream.str()}});
		    charbags.emplace_back(*cs);
		    charbag_map[*cs] = words.size()-1;
		} else
		    words[it->second].emplace_back(line_sstream.str());
	    }
	}
	if (newline == -1)
	    break;
	else
	    line_start = newline+1;
    }

    // Sort the words by charbag hash.

    vector<size_t> hash_sort_order(charbags.size());
    std::iota(hash_sort_order.begin(), hash_sort_order.end(), 0);

    std::sort(hash_sort_order.begin(), hash_sort_order.end(),
	      [&charbags](int a, int b) { return charbags[a].hash() < charbags[b].hash(); });

    apply_permutation_in_place(words, hash_sort_order);
    apply_permutation_in_place(charbags, hash_sort_order);

    //cerr << "Loaded " << count << " dictionary words, " << words.size() << " distinct." << endl;
    return {words, charbags};
}

template <typename Fn>
void forAllAnagrams_iter_last(const vector<CharBag> &dict_charbags,
			      const vector<int> &possible_charbags,
			      const CharBag &charbag, Fn &&f,
			      vector<size_t> &words, size_t start_idx) {
    size_t required_hash = charbag.hash();

    auto compare_charbags = [&dict_charbags, required_hash](const int a, const int b) {
	size_t a_hash, b_hash;

	// FIXME how to do this properly? This is really ugly.
	if (a == -1)
	    a_hash = required_hash;
	else
	    a_hash = dict_charbags[a].hash();
	if (b == -1)
	    b_hash = required_hash;
	else
	    b_hash = dict_charbags[b].hash();

	return a_hash < b_hash;
    };

    {
	auto it = std::lower_bound(possible_charbags.begin() + start_idx,
				   possible_charbags.end(), -1, compare_charbags);

	for (auto ie = possible_charbags.end();
	     it != ie && dict_charbags[*it].hash() == required_hash; ++it)
	    if (charbag == dict_charbags[*it]) {
		words.emplace_back(*it);
		f(words);
		words.pop_back();
		// There's at most one hit, and we found one.
		return;
	    }
    }
}

template <typename Fn>
void forAllAnagrams_iter(const vector<CharBag> &dict_charbags,
			 const vector<int> &old_possible_charbags,
			 const CharBag &charbag, Fn &&f,
			 vector<size_t> &words, size_t start_idx, int curr_len, int max_len) {
    if (curr_len+1 >= max_len) {
	forAllAnagrams_iter_last(dict_charbags, old_possible_charbags, charbag, f,
				 words, start_idx);
	return;
    }
    vector<int> possible_charbags;
    for (int i = start_idx, ie = old_possible_charbags.size(); i<ie; i++)
	if (dict_charbags[old_possible_charbags[i]].isSubsetOf(charbag))
	    possible_charbags.emplace_back(old_possible_charbags[i]);

    for (int i = 0, ie = possible_charbags.size(); i<ie; i++) {
	auto cs = (charbag - dict_charbags[possible_charbags[i]]).value();
	words.emplace_back(possible_charbags[i]);
	if (cs.empty())
	    f(words);
	else
	    forAllAnagrams_iter(dict_charbags, possible_charbags,
				cs, forward<Fn>(f), words, i, curr_len+1, max_len);
	words.pop_back();
    }
}

template <typename Fn>
void forAllAnagrams(const vector<CharBag> &dict_charbags, const CharBag &charbag,
		    int max_len, Fn &&f) {
    vector<size_t> words;
    vector<int> possible_charbags(dict_charbags.size());
    std::iota(possible_charbags.begin(), possible_charbags.end(), 0);
    forAllAnagrams_iter(dict_charbags, possible_charbags,
			charbag, forward<Fn>(f), words, 0, 0, max_len);
}

// The words vector contains vectors of anagram-equivalent words.
// Output all possible combinations of them.
static void outputWords(ostream &stream, const vector<size_t> &word_idxs,
			const vector<vector<string>> &words) {
    int size = word_idxs.size();
    vector<size_t> idxs(size);

    while (true) {
	stream << words[word_idxs[0]][idxs[0]];
	for (int i = 1; i<size; i++)
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
	("dict,d", po::value<string>()->default_value("words.txt"),
	 "dictionary (word list) to use")
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

CharMap generateCharMap(const UnicodeString &input) {
    CharMap charmap;

    unordered_map<UChar, int> char_counts;

    forAllAlpha(input, [&char_counts](UChar c) {
	    ++char_counts[c];
	    return true;
	});

    if (char_counts.size() > MAX_LETTERS) {
	cerr << "Error: More than " << MAX_CHARIDX+1
	     << " different characters in input." << endl;
	exit(1);
    }

    vector<pair<int, UChar>> charmap_with_counts;
    for (const auto &it : char_counts)
	charmap_with_counts.emplace_back(-it.second, it.first);

    std::sort(charmap_with_counts.begin(), charmap_with_counts.end());

    int i = 0;
    for (const auto &p : charmap_with_counts)
	charmap[p.second] = i++;

    return charmap;
}

int main(int argc, char **argv) {
    std::locale::global(std::locale(""));

    auto vm = parse_args(argc, argv);

    UnicodeString input = UnicodeString(vm["sentence"].as<string>().c_str()).toLower();

    CharMap charmap;
    vector<UChar> reverse_charmap;

    charmap = generateCharMap(input);

    auto input_charbag = CharBag::fromUString(input, charmap).value();

    vector<vector<string>> dict_words;
    vector<CharBag> dict_charbags;
    tie(dict_words, dict_charbags) = loadDictionary(
	vm["dict"].as<string>(), input_charbag, charmap);

    forAllAnagrams(dict_charbags, input_charbag, vm["len"].as<int>(),
		   [&dict_words](const vector<size_t> &word_idxs) {
	    assert(!word_idxs.empty());
	    outputWords(cout, word_idxs, dict_words);
	});
}
