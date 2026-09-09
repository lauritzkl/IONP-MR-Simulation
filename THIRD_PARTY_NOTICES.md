# Third-party notices

The project's [MIT license](LICENSE) does not replace the copyright notices
of third-party components. Their original notices remain in the distributed
files.

| Component | Distributed location | Copyright and license |
| --- | --- | --- |
| [CSV for C++ 2.1.3](https://github.com/vincentlaucsb/csv-parser) | `external/csv-parser/csv.hpp` | Copyright (c) 2017–2020 Vincent La; MIT |
| [mio](https://github.com/mandreyel/mio), embedded in the CSV header | `external/csv-parser/csv.hpp` | Copyright 2017 https://github.com/mandreyel; MIT |

The CSV header is retained for source completeness; the current simulator
build does not compile or include it. Both full license notices are preserved
inside that header.

The build uses the compiler's C++ standard library and an installed OpenMP
implementation. Optional plotting uses separately installed NumPy and
Matplotlib. Those dependencies are not bundled in this source release and
retain their respective licenses.
