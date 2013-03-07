#include <cstdio>
#include <cstring>
using namespace std;

const int MAXM = 30, KIND = 26;
int m;

struct node {
    char *s;
    int prefix;   /*作为某些字符串前缀的次数*/
    bool isword;  /*从根到此是一个单词*/
    node *next[KIND]; /*26个英文字母子树，节点太多了，有点浪费哈*/
    node () {
        s = NULL;
        prefix = 0;
        isword = false;
        memset(next, 0, sizeof(next));
    }
} *root;
void insert(node *root, char *s) {
    node *p = root;
    for (int i=0; s[i]; i++) {
        int x = s[i] - 'a';
        p->s = s + i;
        if (p->next[x] == NULL)
            p->next[x] = new node;
        p = p->next[x];
        p->prefix++;
    }
    p->isword = true;
}
bool del(node *root, char *s) {
    node *p = root;
    for (int i=0; s[i]; i++) {
        int x = s[i] - 'a';
        if (p->next[x] == NULL)
            return false;
        p = p->next[x];
    }
    if (p->isword)
        p->isword = false;
    else return false;
    return true;
}
bool search(node *root, char *s) {
    node *p = root;
    for (int i=0; s[i]; i++) {
        int x = s[i] - 'a';
        if (p->next[x] == NULL)
            return false;
        p = p->next[x];
    }
    return p->isword;
}
int count(node *root, char *s) {
    node *p = root;
    for (int i=0; s[i]; i++) {
        int x = s[i] - 'a';
        if (p->next[x] == NULL)
            return 0;
        p = p->next[x];
    }
    return p->prefix;
}
int main() {
    m = 0;
    root = new node;
    char s[MAXM];
    while (gets(s)) {
        if (strcmp(s, "") == 0)
            break;
        insert(root, s);
    }
    while (gets(s)) {
        printf("%d\n", count(root, s));
        if (search(root, s))
            printf("This is a word in Trie.\n");
        else printf("This is not a word in Trie.\n");
    }

    return 0;
}
