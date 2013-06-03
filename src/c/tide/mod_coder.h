// Benjamin Diament

#ifndef MOD_CODER_H
#define MOD_CODER_H

class ModCoder {
 public:
  typedef int Mod;

  void Init(int num_unique_deltas) {
    log_unique_deltas_ = IntLog(num_unique_deltas);
    mask_ = (1 << log_unique_deltas_) - 1;
  }

  int EncodeMod(int aa_index, int unique_delta_index) {
    return (aa_index << log_unique_deltas_) + unique_delta_index;
  }

  void DecodeMod(int code, int* aa_index, int* unique_delta_index) {
    *aa_index = code >> log_unique_deltas_;
    *unique_delta_index = code & mask_;
  }

 private:
  static int IntLog(int x) {
    if (x <= 1)
      return x;
    int res = 0;
    for (--x; x > 0; x >>= 1)
      ++res;
    return res;
  }

  int log_unique_deltas_;
  int mask_;
};

#endif // MOD_CODER_H
