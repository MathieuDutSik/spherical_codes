// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStab.h"
#include "POLY_RecursiveDualDesc.h"
// clang-format on

template <typename T> MyMatrix<T> DropZeroColumn(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  std::vector<int> lnz;
  auto is_nz = [&](int const &iCol) -> bool {
    for (int iRow = 0; iRow < nbRow; iRow++) {
      if (M(iRow, iCol) != 0) {
        return true;
      }
    }
    return false;
  };
  for (int iCol = 0; iCol < nbCol; iCol++) {
    if (is_nz(iCol)) {
      lnz.push_back(iCol);
    }
  }
  MyMatrix<T> Mred(nbRow, lnz.size());
  int pos = 0;
  for (auto &eCol : lnz) {
    for (int iRow = 0; iRow < nbRow; iRow++) {
      Mred(iRow, pos) = M(iRow, eCol);
    }
    pos += 1;
  }
  return Mred;
}

template <typename T, typename Tgroup>
void process_entry_type(std::string const &FileGram,
                        std::string const &FileCode1,
                        std::string const &FileCode2) {
  using TintGroup = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  MyMatrix<T> preCODE1 = ReadMatrixFile<T>(FileCode1);
  MyMatrix<T> preCODE2 = ReadMatrixFile<T>(FileCode2);
  int nbEnt = preCODE1.rows();
  int nbCol = preCODE1.cols();
  std::cerr << "nbEnt=" << nbEnt << " nbCol=" << nbCol << "\n";
  MyMatrix<T> CODE1 = DropZeroColumn(preCODE1);
  MyMatrix<T> CODE2 = DropZeroColumn(preCODE2);
  int dim = CODE1.cols();
  auto get_gram_mat = [&]() -> MyMatrix<T> {
    if (FileGram == "identity") {
      return IdentityMat<T>(dim);
    } else {
      return ReadMatrixFile<T>(FileGram);
    }
  };
  MyMatrix<T> GramMat = get_gram_mat();
  //
  int rnk = RankMat(CODE1);
  std::cerr << "nbEnt=" << nbEnt << " dim=" << dim << " rnk=" << rnk << "\n";
  if (dim != rnk) {
    std::cerr << "We have dim != rnk\n";
    throw TerminalException{1};
  }
  std::optional<std::vector<Tidx>> opt =
    LinPolytope_Isomorphism_GramMat<T, Tidx>(CODE1, GramMat, CODE2, GramMat, std::cerr);
  if (opt) {
    std::cerr << "They are isomorphic\n";
  } else {
    std::cerr << "They are NOT isomorphic\n";
  }

}

template <typename Tgroup>
void process_entry(std::string const &arith,
                   std::string const &FileGram,
                   std::string const &FileCode1,
                   std::string const &FileCode2) {
  if (arith == "rational") {
    using T = mpq_class;
    return process_entry_type<T, Tgroup>(FileGram, FileCode1, FileCode2);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return process_entry_type<T, Tgroup>(FileGram, FileCode1, FileCode2);
  }
  if (arith == "Qsqrt3") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 3>;
    return process_entry_type<T, Tgroup>(FileGram, FileCode1, FileCode2);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return process_entry_type<T, Tgroup>(FileGram, FileCode1, FileCode2);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 5) {
      std::cerr << "CODE_Analysis [arith] [FileGRAM] [FileCODE1] [FileCODE2]\n";
      std::cerr << "   ---------------\n";
      std::cerr << "arith: rational, Qsqrt2, Qsqrt3, Qsqrt5\n";
      throw TerminalException{1};
    }
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string FileGram = argv[2];
    std::string FileCode1 = argv[3];
    std::string FileCode2 = argv[4];
    process_entry<Tgroup>(arith, FileGram, FileCode1, FileCode2);
    std::cerr << "Normal termination of the program time=" << time << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CODE_Analysis time=" << time << "\n";
    exit(e.eVal);
  }
}
