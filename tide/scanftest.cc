#include<stdio.h>
#include<string.h>
#include<limits.h>
#include<assert.h>

void Report(unsigned int limit, const char* aa, int aa_len,
	    double delta) {
  char buf[aa_len + 1];
  strncpy(buf, aa, aa_len);
  buf[aa_len] = '\0';
  fprintf(stderr, "Limit: %u, Amino acids: '%s', delta: %g\n", limit, buf, delta);
}

#define Error(pos, msg) { \
    fprintf(stderr, "Parsing error at position %d: %s\n", (pos), (msg)); \
    fprintf(stderr, "%s\n", text); \
    for (int i = 0; i < (pos); ++i) fprintf(stderr, " "); \
    fprintf(stderr, "^\n"); \
    return 1; \
  }

int main(int argc, char* argv[]) {
  char* text = argv[1];
  printf("Parsing text: '%s'\n", text);

  int pos = 0;
  while (true) {
    char c;
    int next_pos = -1;
    sscanf(text + pos, " %c%n", &c, &next_pos);
    if (next_pos == -1)
      Error(pos, "Expected modification specification.");
    pos += next_pos - 1;
    unsigned int limit = 0;
    if (c >= '1' && c <= '9') {
      sscanf(text + pos, "%u %n", &limit, &next_pos);
      if (limit == UINT_MAX)
	Error(pos, "Limit too big.");
      pos += next_pos;
    }
    int aa_len = -1, plus_pos = -1, delta_pos = -1, end_pos = -1;
    sscanf(text + pos, "%*[ACDEFGHIKLMNPQRSTVWY]%n %n+ %n%*[0-9.] %n",
	   &aa_len, &plus_pos, &delta_pos, &end_pos);
    if (aa_len == -1)
      Error(pos, "Expected amino acid symbol.");
    assert(plus_pos != -1);
    if (delta_pos == -1)
      Error(pos + plus_pos, "Expected '+' and modification amount.");
    if (end_pos == -1)
      Error(pos + delta_pos, "Expected modification amount.");
    if ((limit == 0) && (aa_len != 1))
      Error(pos, "Static modifications must be specified for one amino acid "
	    "at a time.");
    int confirm_end_pos = -1;
    double delta;
    sscanf(text + pos + delta_pos, "%lg %n", &delta, &confirm_end_pos);
    if (delta_pos + confirm_end_pos != end_pos)
      Error(pos + delta_pos, "Cannot parse modification amount.");
    if (delta <= 0)
      Error(pos + delta_pos, "Modification amount must be positive.");
    Report(limit, text + pos, aa_len, delta);
    pos += end_pos;
    if (text[pos] == '\0')
      break;
    if (text[pos] == ',')
      ++pos;
  }

  return 0;
}
