#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
#include<map>
#include<unordered_set>
#include<unordered_map>
#include<map>
#include<set>
#include<random>
#include<algorithm>
#include<queue>
#include<stack>
#include<deque>
#include<bitset>
#include<cstdio>
#include<cassert>
#include<sstream>
#include<set>

#define int long long int

using namespace std;

class Rational {
private:
    
    int a, b;
    
public:
    
    Rational() : a(0), b(1) {}
    
    Rational(int _a, int _b = 1) : a(_a), b(_b) {
        int g = gcd(abs(a), abs(b));
        a /= g; b /= g;
        if (b < 0) a *= -1;
        b = abs(b);
    }
    
    Rational operator- () const {
        return Rational(-a, b);
    }
    
    Rational operator+ (const Rational& other) const {
        int na = a * other.b + other.a * b;
        int nb = b * other.b;
        int g = gcd(na, nb);
        na /= g; nb /= g;
        return Rational(na, nb);
    }
    
    Rational operator- (const Rational& other) const {
        return *this + (-other);
    }
    
    Rational operator* (const Rational& other) const {
        int na = a * other.a;
        int nb = b * other.b;
        int g = gcd(na, nb);
        na /= g; nb /= g;
        return Rational(na, nb);
    }
    
    Rational operator/ (const Rational& other) const {
        int na = a * other.b;
        int nb = b * other.a;
        int g = gcd(na, nb);
        na /= g; nb /= g;
        return Rational(na, nb);
    }
    
    void operator+= (const Rational& other) {
        *this = *this + other;
    }
    
    void operator-= (const Rational& other) {
        *this = *this - other;
    }
    
    void operator*= (const Rational& other) {
        *this = *this * other;
    }
    
    void operator/= (const Rational& other) {
        *this = *this / other;
    }
    
    string to_string() const {
        if (b == 1) return std::to_string(a);
        else return std::to_string(a) + "/" = std::to_string(b);
    }
    
    friend ostream& operator<< (ostream& out, const Rational& a) {
        return out << a.to_string();
    }
    
    bool operator== (const Rational& other) const {
        return a == other.a && b == other.b;
    }
    
    bool operator!= (const Rational& other) const {
        return !(*this == other);
    }
    
};

class Matrix {
private:
    
    vector<vector<Rational>> table;
    
    pair<int, int> getFirst(int st_row, int st_col) const {
        for (int j = st_col; j < table.front().size(); ++j) {
            for (int i = st_row; i < table.size(); ++i) {
                if (table[i][j] != 0) return {i, j};
            }
        }
        return {table.size(), table.front().size()};
    }
    
    void multiply(int row_i, Rational x) {
        for (int w = 0; w < table.front().size(); ++w) table[row_i][w] *= x;
    }
    
    void add(int row_j, int row_i, Rational x) {
        for (int w = 0; w < table.front().size(); ++w) table[row_j][w] += table[row_i][w] * x;
    }
    
public:
    
    Matrix() {}
    
    Matrix(const vector<vector<Rational>>& a) : table(a) {}
    
    void oneStepGauss(int st_row, int st_col) {
        auto [i, j] = getFirst(st_row, st_col);
        swap(table[st_row], table[i]);
        for (int w = st_row + 1; w < table.size(); ++w) {
            multiply(w, table[st_row][j]);
            add(w, st_row, -table[w][j] / table[st_row][j]);
        }
    }
    
    void getTriangularView() {
        auto p = getFirst(0, 0);
        int now_row = p.first, now_col = p.second;
        while (now_row < table.size() && now_col < table.front().size()) {
            oneStepGauss(now_row, now_col);
            p = getFirst(now_row + 1, now_col + 1);
            now_row = p.first; now_col = p.second;
        }
    }
    
    void getTriangularView_Tex(const string& key) {
        auto p = getFirst(0, 0);
        int now_row = p.first, now_col = p.second;
        while (now_row < table.size() && now_col < table.front().size()) {
            oneStepGauss(now_row, now_col);
            cout << getTex(key);
            p = getFirst(now_row + 1, now_col + 1);
            if (p.first < table.size()) cout << "\\xrightarrow[]{}\n";
            now_row = p.first; now_col = p.second;
        }
    }
    
    void getImprovedSteppedView() {
        auto p = getFirst(0, 0);
        int now_row = p.first, now_col = p.second;
        vector<pair<int, int>> roots;
        while (now_row < table.size() && now_col < table.front().size()) {
            roots.push_back(p);
            oneStepGauss(now_row, now_col);
            p = getFirst(now_row + 1, now_col + 1);
            now_row = p.first; now_col = p.second;
        }
        reverse(roots.begin(), roots.end());
        for (const auto& [i, j] : roots) {
            multiply(i, Rational(1) / table[i][j]);
            for (int w = i - 1; w >= 0; --w) {
                add(w, i, -table[w][j]);
            }
        }
    }
    
    void getImprovedSteppedView_Tex(const string& key) {
        auto p = getFirst(0, 0);
        int now_row = p.first, now_col = p.second;
        vector<pair<int, int>> roots;
        string ans = "";
        ans += getTex(key) + "\\xrightarrow[]{}\n";
        while (now_row < table.size() && now_col < table.front().size()) {
            roots.push_back(p);
            oneStepGauss(now_row, now_col);
            ans += getTex(key) + "\\xrightarrow[]{}\n";
            p = getFirst(now_row + 1, now_col + 1);
            now_row = p.first; now_col = p.second;
        }
        reverse(roots.begin(), roots.end());
        for (const auto& [i, j] : roots) {
            multiply(i, Rational(1) / table[i][j]);
            for (int w = i - 1; w >= 0; --w) {
                add(w, i, -table[w][j]);
            }
            ans += getTex(key) + "\\xrightarrow[]{}\n";
        }
        ans.pop_back(); ans.pop_back(); ans.pop_back(); ans.pop_back();
        while (ans.back() != ']') ans.pop_back();
        cout << ans << "\n";
    }
    
    string getTex(const string& key) const {
        string s = "\\left[\\begin{array}{" + key + "}\n";
        for (const auto & i : table) {
            s += "\t";
            for (int j = 0; j < i.size(); ++j) {
                s += i[j].to_string() + " ";
                if (j < (int)i.size() - 1) s += "& ";
            }
            s += "\\\\ \n";
        }
        s += "\\end{array}\\right]\n";
        return s;
    }
    
    void printTex(const string& key) const {
        cout << getTex(key) << "\n";
    }
    
    void transpose() {
        vector<vector<Rational>> nxt(table.front().size(), vector<Rational>(table.size()));
        for (int i = 0; i < table.front().size(); ++i) {
            for (int j = 0; j < table.size(); ++j) nxt[i][j] = table[j][i];
        }
        table.swap(nxt);
    }
    
    friend istream& operator>> (istream& in, Matrix& a) {
        string line;
        while (getline(in, line)) {
            stringstream ss(line);
            int x;
            vector<Rational> ansLine;
            while (ss >> x) ansLine.push_back(x);
            a.table.push_back(ansLine);
        }
        return in;
    }

    friend ostream& operator<< (ostream& out, Matrix& a) {
        for (const auto& i : a.table) {
            for (auto j : i) out << j << " ";
            out << "\n";
        }
        return out;
    }
};

signed main() {
    Matrix a; cin >> a;
    cout << a << "\n\n";
    a.getImprovedSteppedView_Tex("ccccc|ccc");
    a.transpose();
    cout << a << "\n";
    return 0;
}
