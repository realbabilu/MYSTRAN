# 2026-10-03 — T6 recovery, OP2 and GPSTRESS bundle

Paket ini merangkum selisih dari clone GitHub awal sampai kondisi lokal yang disiapkan untuk commit. Tanggal memakai Asia/Jakarta.

| Item | Nilai |
|---|---|
| Repository | `https://github.com/realbabilu/MYSTRAN.git` |
| Branch | `v18.00.a` |
| Baseline / HEAD saat packaging | `931dde4df56d9c39b43d884bda941f7fd27c6036` |
| File tracked berubah | 38 Fortran files |
| Snapshot Python | 37 files, termasuk referensi dan dependensi impor |
| Salinan file sumber / pendukung / bukti | 136; rincian dalam manifest |
| Commit / push | Belum dilakukan |

## Isi paket

- [CHANGES.md](CHANGES.md): ringkasan lengkap, provenance, keputusan bundle, dan batas validasi.
- [CHANGED_FILES.md](CHANGED_FILES.md): seluruh 38 file berubah dan penjelasan per file, ditambah perubahan Python dengan bukti backup lokal.
- [VALIDATION.md](VALIDATION.md): hasil pengujian dan petunjuk menjalankan ulang.
- [COMMIT_MESSAGE.md](COMMIT_MESSAGE.md): judul dan body commit yang disarankan.
- `files/Source/`: salinan persis versi akhir semua 38 file source yang berubah.
- `SOURCE_FROM_931dde4.patch`: diff tracked source dari baseline GitHub; tidak perlu diaplikasikan ke working tree saat ini karena perubahannya sudah ada.
- `PYTHON_LOCAL_CHANGES/`: diff Python terhadap backup lokal, bukan terhadap GitHub.
- `validation/python/` dan `validation/test/`: benchmark, referensi Python final, helper, skrip validasi, serta deck.
- `evidence/`: hasil JSON, laporan perbandingan, catatan resume, dan log build/run yang dipilih.
- `history/one_time_migrations/`: catatan skrip migrasi lama; jangan dijalankan ulang sebagai pengujian.
- `FILE_MANIFEST.csv` / `manifest.json`: asal setiap salinan, jenis/provenance, ukuran dan SHA-256.
- `CHECKSUMS.sha256`: integritas seluruh file bundle kecuali file checksum itu sendiri.
- `tools/verify_bundle.py`: verifikasi checksum dan kesesuaian snapshot source dengan working tree.
- `tools/restore_validation.py`: menyalin aset validasi kembali ke layout sibling `python`/`test` bila dibutuhkan pada clone lain.

Satu bundle digunakan karena formula T6, ukuran array, pemulihan stress/force, GPSTRESS dan OP2 berbagi kontrak data. Memisahkan commit tanpa pengujian ulang pada setiap tahap dapat menghasilkan kombinasi source yang tidak konsisten. Di dalam dokumentasi, perubahan tetap dibagi menurut fungsi agar review mudah.

`Source/` di repo adalah kode produksi. Folder `files/` tidak menjadi input tambahan CMake. Executable, object/module compiler, cache Python, DLL, PAT/kredensial, serta file OP2/F06 besar tidak disertakan. File tersebut dapat dibentuk ulang melalui skrip validasi.
