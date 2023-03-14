EXECUTABLES_REQUIREMENTS = {"meryl": {"executable": "meryl",
                                      "user_path": "MERYL_PATH"},
                            "gffread": {"executable": "gffread",
                                        "user_path": "GFFREAD_PATH"},
                            "samtools": {"executable": "samtools",
                                        "user_path": "SAMTOOLS_PATH"},
                            "bedtools": {"executable": "bedtools",
                                         "user_path": "BEDTOOLS_PATH"}}

MAGIC_NUMS_COMPRESSED = [b'\x1f\x8b\x08', b'\x42\x5a\x68', 
                         b'\x50\x4b\x03\x04']