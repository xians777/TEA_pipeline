//============================================================================
// Name        : TEA.hpp
// Author      : Kyu Park
// Version     :
// Copyright   : Lee Lab @ Boston Children's Hospital
// Description : TEA 2.0
//============================================================================


#ifndef TEA_TEA_HPP_
#define TEA_TEA_HPP_

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/range/adaptor/reversed.hpp>

#include <api/BamAux.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "../castle/ParallelRunner.hpp"
#include "../castle/StringUtils.hpp"
#include "../castle/IOUtils.hpp"
#include "../castle/MathUtils.hpp"
#include "../meerkat/BWACaller.hpp"
#include "../meerkat/BlockBoundary.hpp"
#include "../meerkat/preprocess/ReadGroup.h"
//#include "../meerkat/dre/ParallelDiscordExtractor.hpp"
#include "../meerkat/BlockBoundary.hpp"
#include "TEAOptionParser.hpp"
#include "RAMEntry.hpp"


namespace tea {
using namespace BamTools;


class TEA {
	struct Read {
		struct Clip {
			string	name;
			uint16_t mapq;
			char	cigar;
			string	bases;
			string	quals;

			Clip() : name(""), mapq(0), cigar(' '), bases(""), quals("") {
			}

			Clip(const Clip& other) {
				name = other.name;
				mapq = other.mapq;
				cigar = other.cigar;
				bases = other.bases;
				quals = other.quals;
			}

			Clip& operator=(const Clip& other) {
				name = other.name;
				mapq = other.mapq;
				cigar = other.cigar;
				bases = other.bases;
				quals = other.quals;

				return *this;
			}

			string fastq(void) {
				return (boost::format("@%s\n%s\n+\n%s\n") % name % bases % quals).str();
			}

			bool empty(void) {
				if (name == "" || bases == "" || quals == "")
					return true;
				else
					return false;
			}
		};

		string	name;
		uint16_t mapq;
		string	bases;
		string	quals;
		Clip 	head;
		Clip 	tail;
		Clip	fullmap;
		bool	is_first;
		bool	is_second;
		bool	is_primary;
		bool	is_supplementary;

		Read() : name(), mapq(0), bases(), quals(),
				head(), tail(), fullmap(),
				is_first(false), is_second(false), is_primary(false), is_supplementary(false)
				{
		}

		Read(const Read& other) {
			name = other.name;
			mapq = other.mapq;
			bases = other.bases;
			quals = other.quals;
			head = other.head;
			tail = other.tail;
			fullmap = other.fullmap;
			is_first = other.is_first;
			is_second = other.is_second;
			is_primary = other.is_primary;
			is_supplementary = other.is_supplementary;
		}

		Read(const BamAlignment& al) {

			bases = al.QueryBases;
			quals = al.Qualities;

			if (bases.size() == 0 || quals.size() == 0 || bases.size() != quals.size()) {
				return;
			}

			name = al.Name;
			mapq = al.MapQuality;
			is_first = al.IsFirstMate();
			is_second = al.IsSecondMate();
			is_primary = al.IsPrimaryAlignment();
			is_supplementary = al.IsSupplementary();

			string debug_name = "HWI-ST1221:217:D1R5UACXX:1:2105:9090:32286";
			bool debug = (name == debug_name);

			auto& cigars = al.CigarData;
			auto cigars_front_type 		= cigars.front().Type;
			auto cigars_front_length 	= cigars.front().Length;
			auto cigars_back_type 		= cigars.back().Type;
			auto cigars_back_length 	= cigars.back().Length;

			if ('M' == cigars_front_type && 'M' == cigars_back_type ) {
				fullmap.name = name;
				fullmap.mapq = mapq;
				fullmap.cigar = 'M';
				fullmap.bases = bases;
				fullmap.quals = quals;
			}
			else {
				if (al.IsReverseStrand() && cigars.size() >= 2) {
					if ( ('S' == cigars_front_type || 'H' == cigars_front_type)
							&& ('S' == cigars_back_type || 'H' == cigars_back_type) ) {

						if (cigars_front_length > cigars_back_length) {
							if ('H' == cigars_front_type) cigars_front_length = 0;
							if ('H' == cigars_back_type) cigars_back_length = 0;

							if (debug) {
								cout << "[Read(BamAlignment)] " << name << "\t" << cigars_front_type << "\treverse, clips at both side, keep front clip \n";
							}

							head.name = tail.name = name;
							head.mapq = tail.mapq = mapq;
							head.cigar = cigars_front_type;
							tail.cigar = 'M';


							head.bases = bases.substr(0, cigars_front_length);
							head.quals = quals.substr(0, cigars_front_length);


							tail.bases = bases.substr(cigars_front_length, al.QueryBases.size() - cigars_front_length - cigars_back_length);
							tail.quals = quals.substr(cigars_front_length, al.Qualities.size() - cigars_front_length - cigars_back_length);

						}

						else if (cigars_front_length < cigars_back_length) {
							if ('H' == cigars_front_type) cigars_front_length = 0;
							if ('H' == cigars_back_type) cigars_back_length = 0;

							if (debug) {
								cout << "[Read(BamAlignment)] " << name << "\t" << cigars_back_type << "\treverse, clips at both side, keep back clip \n";
							}

							head.name = tail.name = name;
							head.mapq = tail.mapq = mapq;

							head.cigar = 'M';
							head.bases = bases.substr(cigars_front_length, al.QueryBases.size() - cigars_front_length - cigars_back_length);
							head.quals = quals.substr(cigars_front_length, al.Qualities.size() - cigars_front_length - cigars_back_length);

							tail.cigar = cigars_back_type;
							tail.bases = bases.substr(al.QueryBases.size() - cigars_back_length);
							tail.quals = quals.substr(al.Qualities.size() - cigars_back_length);
						}
						else {
							cout << "same size clippings at both sides: " << str() << "\n";
						}
					}

					else if ('S' == cigars_front_type || 'H' == cigars_front_type) {
						if ('H' == cigars_front_type) cigars_front_length = 0;
						if ('H' == cigars_back_type) cigars_back_length = 0;

						if (debug) {
							cout << "[Read(BamAlignment)] " << name << "\t" << cigars_front_type << "\treverse, single front clip \n";
						}

						head.name = tail.name = name;
						head.mapq = tail.mapq = mapq;

						head.cigar = cigars_front_type;
						head.bases = bases.substr(0, cigars_front_length);
						head.quals = quals.substr(0, cigars_front_length);

						tail.cigar = 'M';
						tail.bases = bases.substr(cigars_front_length);
						tail.quals = quals.substr(cigars_front_length);
					}

					else if ('S' == cigars_back_type || 'H' == cigars_back_type) {
						if ('H' == cigars_front_type) cigars_front_length = 0;
						if ('H' == cigars_back_type) cigars_back_length = 0;

						if (debug) {
							cout << "[Read(BamAlignment)] " << name << "\t" << cigars_back_type << "\treverse, single back clip \n";
						}

						head.name = tail.name = name;
						head.mapq = tail.mapq = mapq;

						head.cigar = 'M';
						head.bases = bases.substr(0, al.QueryBases.size() - cigars_back_length);
						head.quals = quals.substr(0, al.Qualities.size() - cigars_back_length);

						tail.cigar = cigars_back_type;
						tail.bases = bases.substr(al.QueryBases.size() - cigars_back_length);
						tail.quals = quals.substr(al.Qualities.size() - cigars_back_length);
					}
				}
				else if (!al.IsReverseStrand() && cigars.size() >= 2) {
					if ( ('S' == cigars_front_type || 'H' == cigars_front_type)
							&& ('S' == cigars_back_type || 'H' == cigars_back_type) ) {

						if (cigars_front_length > cigars_back_length) {
							if ('H' == cigars_front_type) cigars_front_length = 0;
							if ('H' == cigars_back_type) cigars_back_length = 0;

							if (debug) {
								cout << "[Read(BamAlignment)] " << name << "\t" << cigars_front_type << "\tforward, clips at both side, keep front clip \n";
							}

							head.name = tail.name = name;
							head.mapq = tail.mapq = mapq;

							head.cigar = 'M';
							head.bases = bases.substr(cigars_front_length, al.QueryBases.size() - cigars_front_length - cigars_back_length);
							head.quals = quals.substr(cigars_front_length, al.QueryBases.size() - cigars_front_length - cigars_back_length);

							tail.cigar = cigars_front_type;
							tail.bases = bases.substr(0, cigars_front_length);
							tail.quals = quals.substr(0, cigars_front_length);
						}

						else if (cigars_front_length < cigars_back_length) {
							if ('H' == cigars_front_type) cigars_front_length = 0;
							if ('H' == cigars_back_type) cigars_back_length = 0;

							if (debug) {
								cout << "[Read(BamAlignment)] " << name << "\t" << cigars_back_type << "\tforward, clips at both side, keep back clip \n";
							}

							head.name = tail.name = name;
							head.mapq = tail.mapq = mapq;

							head.cigar = cigars_back_type;
							head.bases = bases.substr(al.QueryBases.size() - cigars_back_length);
							head.quals = quals.substr(al.Qualities.size() - cigars_back_length);

							tail.cigar = 'M';
							tail.bases = bases.substr(cigars_front_length, al.QueryBases.size() - cigars_front_length - cigars_back_length);
							tail.quals = quals.substr(cigars_front_length, al.Qualities.size() - cigars_front_length - cigars_back_length);
						}

						else {
							cout << "same size clippings at both sides: " << str() << "\n";
						}
					}

					else if ('S' == cigars_front_type || 'H' == cigars_front_type) {
						if ('H' == cigars_front_type) cigars_front_length = 0;
						if ('H' == cigars_back_type) cigars_back_length = 0;

						if (debug) {
							cout << "[Read(BamAlignment)] " << name << "\t" << cigars_front_type << "\tforward, single front clip \n";
						}

						head.name = tail.name = name;
						head.mapq = tail.mapq = mapq;

						head.cigar = 'M';
						head.bases = bases.substr(cigars_front_length);
						head.quals = quals.substr(cigars_front_length);

						tail.cigar = cigars_front_type;
						tail.bases = bases.substr(0, cigars_front_length);
						tail.quals = quals.substr(0, cigars_front_length);
					}

					else if ('S' == cigars_back_type || 'H' == cigars_back_type) {
						if ('H' == cigars_front_type) cigars_front_length = 0;
						if ('H' == cigars_back_type) cigars_back_length = 0;

						if (debug) {
							cout << "[Read(BamAlignment)] " << name << "\t" << cigars_back_type << "\tforward, single back clip \n";
						}

						head.name = tail.name = name;
						head.mapq = tail.mapq = mapq;

						head.cigar = cigars_back_type;
						head.bases = bases.substr(al.QueryBases.size() - cigars_back_length);
						head.quals = quals.substr(al.Qualities.size() - cigars_back_length);

						tail.cigar = 'M';
						tail.bases = bases.substr(0, al.QueryBases.size() - cigars_back_length);
						tail.quals = quals.substr(0, al.Qualities.size() - cigars_back_length);
					}
				}
				else {
					cout << "[Read(BamAlignment)] unknown case\t" << name << "\t" << bases << "\n";
				}
			}

			if (debug) {
				cout << "[Read(BamAlignment)] " << name << "\t" << head.bases << "\t" << tail.bases << "\n";
			}
		}

		string str(void) {
			return (boost::format("%s\tbases:%s\tfullmap:%s %s\thead:%s %s\ttail:%s %s") % name % bases % fullmap.cigar % fullmap.bases % head.cigar % head.bases % tail.cigar % tail.bases ).str();
		}

		bool m_clip_as_fullmap(void) {
			if (!is_fullmap()) {
				if (head.cigar == 'M' && tail.cigar != 'M') {
					fullmap = head;
					head = Clip();
					tail = Clip();
					return true;
				}
				else if (tail.cigar == 'M' && head.cigar != 'M') {
					fullmap = tail;
					head = Clip();
					tail = Clip();
					return true;
				}
				else return false;
			}
			else return false;
		}

		bool is_fullmap(void) {
			if (!fullmap.empty() && head.empty() && tail.empty())
				return true;
			else
				return false;
		}
		bool empty(void) {
			if ((bases == "" || quals == "") && fullmap.empty() && head.empty() && tail.empty())
				return true;
			else
				return false;
		}
	};

	struct ReadPair {
		string name;
		pair<Read, Read> first_mate;
		pair<Read, Read> second_mate;
		Read first_read;
		Read second_read;
		Read::Clip first_clip;
		Read::Clip second_clip;

		ReadPair() : name(), first_mate(), second_mate(),
				first_read(), second_read(), first_clip(), second_clip() {
		}

		const uint16_t LOW_MQ = 10;

		string debug_name = "HWI-ST1221:217:D1R5UACXX:1:2105:9090:32286";

		void addRead(Read& r) {
			bool debug = (r.name == debug_name);;
			name = r.name;

			if (r.is_first && !r.is_second ) {
				if (r.is_supplementary || !r.is_primary) {
					if (!first_mate.second.empty()) {
						cout << "[addRead(Read)] overwriting in first_mate.second " << first_mate.second.is_primary << " " << first_mate.second.is_supplementary << " " << first_mate.second.str() << "\n";
						debug = false;
					}
					first_mate.second = r;
					if (debug) {
						cout << "[addRead(Read)] first_mate.first: " << " " << first_mate.first.is_primary << " " << first_mate.first.is_supplementary  << " " << first_mate.first.str() << " \n";
						cout << "[addRead(Read)] first_mate.second: " << " " << first_mate.second.is_primary << " " << first_mate.second.is_supplementary  << " " << first_mate.second.str() << " \n";
						cout << "[addRead(Read)] second_mate.first: " << " " << second_mate.first.is_primary << " " << second_mate.first.is_supplementary  << " " << second_mate.first.str() << " \n";
						cout << "[addRead(Read)] second_mate.second: " << " " << second_mate.second.is_primary << " " << second_mate.second.is_supplementary  << " " << second_mate.second.str() << " \n";
					}
				}
				else {
					if (!first_mate.first.empty()) {
						cout << "[addRead(Read)] overwriting in first_mate.first " << first_mate.first.is_primary << " " << first_mate.first.is_supplementary << " " << first_mate.first.str() << "\n";
						debug = false;
					}
					first_mate.first = r;
					if (debug) {
						cout << "[addRead(Read)] first_mate.first: " << " " << first_mate.first.is_primary << " " << first_mate.first.is_supplementary  << " " << first_mate.first.str() << " \n";
						cout << "[addRead(Read)] first_mate.second: " << " " << first_mate.second.is_primary << " " << first_mate.second.is_supplementary  << " " << first_mate.second.str() << " \n";
						cout << "[addRead(Read)] second_mate.first: " << " " << second_mate.first.is_primary << " " << second_mate.first.is_supplementary  << " " << second_mate.first.str() << " \n";
						cout << "[addRead(Read)] second_mate.second: " << " " << second_mate.second.is_primary << " " << second_mate.second.is_supplementary  << " " << second_mate.second.str() << " \n";
					}
				}
			}

			else if (!r.is_first && r.is_second ) {
				if (r.is_supplementary || !r.is_primary) {
					if (!second_mate.second.empty()) {
						cout << "[addRead(Read)] overwriting in second_mate.second " << second_mate.second.is_primary << " " << second_mate.second.is_supplementary << " " << second_mate.second.str() << "\n";
						debug = false;
					}
					second_mate.second = r;
					if (debug) {
						cout << "[addRead(Read)] first_mate.first: " << " " << first_mate.first.is_primary << " " << first_mate.first.is_supplementary  << " " << first_mate.first.str() << " \n";
						cout << "[addRead(Read)] first_mate.second: " << " " << first_mate.second.is_primary << " " << first_mate.second.is_supplementary  << " " << first_mate.second.str() << " \n";
						cout << "[addRead(Read)] second_mate.first: " << " " << second_mate.first.is_primary << " " << second_mate.first.is_supplementary  << " " << second_mate.first.str() << " \n";
						cout << "[addRead(Read)] second_mate.second: " << " " << second_mate.second.is_primary << " " << second_mate.second.is_supplementary  << " " << second_mate.second.str() << " \n";
					}
				}
				else {
					if (!second_mate.first.empty()) {
						cout << "[addRead(Read)] overwriting in second_mate.first " << second_mate.first.is_primary << " " << second_mate.first.is_supplementary << " " << second_mate.first.str() << "\n";
						debug = false;
					}
					second_mate.first = r;
					if (debug) {
						cout << "[addRead(Read)] first_mate.first: " << " " << first_mate.first.is_primary << " " << first_mate.first.is_supplementary  << " " << first_mate.first.str() << " \n";
						cout << "[addRead(Read)] first_mate.second: " << " " << first_mate.second.is_primary << " " << first_mate.second.is_supplementary  << " " << first_mate.second.str() << " \n";
						cout << "[addRead(Read)] second_mate.first: " << " " << second_mate.first.is_primary << " " << second_mate.first.is_supplementary  << " " << second_mate.first.str() << " \n";
						cout << "[addRead(Read)] second_mate.second: " << " " << second_mate.second.is_primary << " " << second_mate.second.is_supplementary  << " " << second_mate.second.str() << " \n";
					}
				}
			}
			else {
				cout << "addRead: no where to add \n";
			}
		}

		void setClips() {
			string debug_name = "HWI-ST1221:217:D1R5UACXX:1:2105:9090:32286";
			bool debug = (name == debug_name);
			if (debug) cout << "[ReadPair.setClips] " << name << "\n";

			setFirstAndSecondReads();
			setFirstAndSecondClips();
		}

		void setFirstAndSecondReads() {
			string debug_name = "HWI-ST1221:217:D1R5UACXX:1:2105:9090:32286";
			bool debug = (name == debug_name);

			//if (!first_mate.first.empty() && first_mate.second.empty()) {
			if (!first_mate.first.empty() && first_mate.second.empty()) {
				if (first_mate.first.is_fullmap()) {
					first_read = first_mate.first;
					if (debug) cout << "[ReadPair.setFirstAndSecondReads] first_mate.first:"
							<< first_read.str() << "\n";
				}
				else if (first_mate.first.m_clip_as_fullmap()) {
					first_read = first_mate.first;
					if (debug) cout << "[ReadPair.setFirstAndSecondReads] m_clip_as_fullmap first_mate.first:"
							<< first_read.str() << "\n";
				}
			}
			else if (first_mate.first.empty() && !first_mate.second.empty()) {
				first_read = first_mate.second;
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] first_mate.second:" << first_read.str() << "\n";
			}
			else if (!first_mate.first.empty() && !first_mate.second.empty()) {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] first_mate has two reads \n";
			}
			else if (first_mate.first.empty() && first_mate.second.empty()) {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] first_mate both empty \n";
			}
			else {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] first_mate \n";
			}

			if (!second_mate.first.empty() && second_mate.second.empty()) {
				if (second_mate.first.is_fullmap()) {
					second_read = second_mate.first;
					if (debug) cout << "[ReadPair.setFirstAndSecondReads] second_mate.first:"
							<< second_read.str() << "\n";
				}
				else if (second_mate.first.m_clip_as_fullmap()) {
					second_read = second_mate.first;
					if (debug) cout << "[ReadPair.setFirstAndSecondReads] m_clip_as_fullmap second_mate.first:"
							<< second_read.str() << "\n";
				}
			}
			else if (second_mate.first.empty() && !second_mate.second.empty()) {
				second_read = second_mate.second;
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] second_mate.second:" << second_read.str() << "\n";
			}
			else if (!second_mate.first.empty() && !second_mate.second.empty()) {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] second_mate has two reads \n";
			}
			else if (second_mate.first.empty() && second_mate.second.empty()) {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] second_mate both empty \n";
			}
			else {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] second_mate \n";
			}

			if (first_mate.first.is_fullmap() && first_mate.second.is_fullmap()) {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] Two fullmaps at first mate:" << first_mate.first.name << "\n";
			}
			else if (second_mate.first.is_fullmap() && second_mate.second.is_fullmap()) {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] Two fullmaps at second mate:" << second_mate.first.name << "\n";
			}
			else if (first_mate.first.is_fullmap() && !first_mate.second.empty()) {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] Fullmap with secondary:" << first_mate.first.name << "\n";
			}
			else if (second_mate.first.is_fullmap() && !second_mate.second.empty()) {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] Fullmap with secondary:" << second_mate.first.name << "\n";
			}
			else {
				if (debug) cout << "[ReadPair.setFirstAndSecondReads] both_mates \n";
			}
		}

		void setFirstAndSecondClips() {
			string debug_name = "HWI-ST1221:217:D1R5UACXX:1:2105:9090:32286";

			bool debug = false;
			debug = (name == debug_name);

			if (!first_read.empty() && !second_read.empty()) {
				if (first_read.is_fullmap() && !second_read.is_fullmap()) {
					setClipsWithFullmapAndNonfullmapReads(first_read, second_read, true);
					if (debug) cout << "[setFirstAndSecondClips] setClipsWithFullmapAndNonfullmapReads\n";
				}
				else if (second_read.is_fullmap() && !first_read.is_fullmap()) {
					setClipsWithFullmapAndNonfullmapReads(second_read, first_read, false);
					if (debug) cout << "[setFirstAndSecondClips] setClipsWithFullmapAndNonfullmapReads\n";
				}
				else if (!first_read.is_fullmap() && !second_read.is_fullmap()) {
					cout << "[setFirstAndSecondClips] No fullmap reads: " << name << "\n";
					return;
				}
				else if (first_read.is_fullmap() && second_read.is_fullmap()) {
					first_clip = first_read.fullmap;
					second_clip = second_read.fullmap;
					if (debug) cout << "[setFirstAndSecondClips] Fullmap reads at both sides: " << name << "\n";
					return;
				}
				else {
					cout << "[setFirstAndSecondClips] unknown case: " << name << "\n";
					return;
				}
			}
			else if (first_read.is_fullmap() && second_read.empty() && second_mate_has_two_reads()) {
				setClipsWithFullmapReadAndTwoReadMate(first_read, second_mate, true);
				if (debug) cout << "[setFirstAndSecondClips] setClipsWithFullmapReadAndTwoReadMate\n";
			}
			else if (second_read.is_fullmap() && first_read.empty() && first_mate_has_two_reads()) {
				setClipsWithFullmapReadAndTwoReadMate(second_read, first_mate, false);
				if (debug) cout << "[setFirstAndSecondClips] setClipsWithFullmapReadAndTwoReadMate\n";
			}
			else {
				cout << "[setFirstAndSecondClips] unknown case " << first_read.name << "\t" << second_read.name << "\n";
			}
		}

		void setClipsWithFullmapReadAndTwoReadMate(
				Read read, pair<Read, Read> mate, bool is_fullmap_in_first) {
			Read::Clip head;
			Read::Clip tail;
			Read::Clip mate_clip;

			string debug_name = "HWI-ST1221:217:D1R5UACXX:1:2105:9090:32286";
			bool debug = (name == debug_name);

			if ( mate.first.head.cigar == 'M' && mate.second.tail.cigar == 'M') {
				head = mate.first.head;
				tail = mate.second.tail;
			}
			else if ( mate.second.head.cigar == 'M' && mate.first.tail.cigar == 'M') {
				head = mate.second.head;
				tail = mate.first.tail;
			}
			else {
				cout << "[setClipsWithFullmapReadAndTwoReadMate] unknown case " << read.name << "\n";
				return;
			}

			if (debug) cout << "head " << head.name << " " << head.bases << "\n";
			if (debug) cout << "tail " << tail.name << " " << tail.bases << "\n";

			if (read.mapq == 0) {
				if (head.mapq < tail.mapq) {
					mate_clip = tail;
					if (debug) cout << "[Bc4G] - Return 100M, Tail(M)\n" << tail.fastq() << " \n";
				}
				else {
					if (debug) cout << "[Bc4H] known case to skip" << read.str() << " \n";
					return;
				}
			}
			else if (head.mapq == tail.mapq) {
				mate_clip = tail;
				if (debug) cout << "[Bd6] - Return 100M, Tail(M)\n";
			}
			else {
				if (read.mapq < LOW_MQ) {
					mate_clip = head.mapq > tail.mapq ? head : tail;
					if (debug) cout << "[Be7J] - Return 100M, Higher MQ (M)\n";
				}
				else {
					mate_clip = head.mapq > tail.mapq ? tail : head;
					if (debug) cout << "[BeK] - Return 100M, Lower MQ (M)\n";
				}
			}

			if (is_fullmap_in_first) {
				first_clip = read.fullmap;
				second_clip = mate_clip;
			}
			else {
				first_clip = mate_clip;
				second_clip = read.fullmap;
			}

			if (debug) cout << "first" << first_clip.fastq() << "second" << second_clip.fastq();

		}


		void setClipsWithFullmapAndNonfullmapReads(Read fullmap, Read nonfullmap, bool is_fullmap_in_first) {

			Read::Clip nonfullmap_clip;

			string debug_name = "HWI-ST1221:217:D1R5UACXX:1:2105:9090:32286";
			bool debug = (name == debug_name);

			if (debug) cout << nonfullmap.head.fastq();
			if (debug) cout << nonfullmap.tail.fastq();

			if (fullmap.mapq >= LOW_MQ) {
				if ('M' == nonfullmap.head.cigar ) {
					if (fullmap.mapq < LOW_MQ) {
						nonfullmap_clip = nonfullmap.head;
						if (debug) {
							cout << "[Aa1Cf] " << nonfullmap_clip.bases << "\t" << nonfullmap.str() << " \n";
						}
					}
					else {
						cout << "[Aa1Dg] known case to skip:" << fullmap.name << " \n";
						return;
					}

				}
				else if ('M' == nonfullmap.tail.cigar ) {
					if (fullmap.mapq >= LOW_MQ) {
						nonfullmap_clip = nonfullmap.head;
						if (debug) {
							cout << "[Aa2Ei] " << nonfullmap_clip.bases << "\t" << nonfullmap.str() << " \n";
						}
					}
					else {
						cout << "[Aa2Fj] known case to skip:" << fullmap.name << " \n";
						return;
					}

				}
				else {
					cout << "[setClipsWithFullmapAndNonfullmapReads] [Aa] unknown case" << fullmap.name << "\n";
					return;
				}
			}
			else {
				nonfullmap_clip = nonfullmap.tail;
				if (debug) {
					cout << "[Ab3] " << nonfullmap_clip.bases << "\t" << nonfullmap.str() << " \n";
				}
			}

			if (is_fullmap_in_first) {
				first_clip = fullmap.fullmap;
				second_clip = nonfullmap_clip;
			}
			else {
				first_clip = nonfullmap_clip;
				second_clip = fullmap.fullmap;
			}

			if (debug) cout << nonfullmap.str() << "\n";
			if (debug) cout << first_clip.fastq();
			if (debug) cout << second_clip.fastq();

		}

		bool is_first_mate_empty() {
			if (!first_mate.first.empty() || !first_mate.second.empty())
				return false;
			else
				return true;
		}

		bool is_second_mate_empty() {
			if (!second_mate.first.empty() || !second_mate.second.empty())
				return false;
			else
				return true;
		}

		bool first_mate_has_two_reads() {
			if (first_mate.first.empty() || first_mate.second.empty())
				return false;
			else
				return true;
		}

		bool second_mate_has_two_reads() {
			if (second_mate.first.empty() || second_mate.second.empty())
				return false;
			else
				return true;
		}

		bool has_both_clips() {
			if ( first_clip.empty() || second_clip.empty() )
				return false;
			else
				return true;
		}

	};


public:
	TEA();
	~TEA();

	void set_option_parser(const TEAOptionParser& options);
	void preprocess();
	void preprocess_v();
	void preprocess_u();
	void run_rid();
	void run_vid();
	void run_uid();

	void append_contig();
	void append_contig_alt();

	// transduction
	void run_transduction();
	void create_tea_transduction(const string& out_path, const string& in_path);
	void create_discord_bed(const string& out_path, const string& in_path);
	void read_ram_read_ids(boost::unordered_set<string>& ram_read_ids);
	void read_ram_ids(boost::unordered_set<string>& ram_ids, const string& in_path, boost::unordered_map<string, int32_t>& cluster_entries);
	void create_umm_tmp_from_cluster_raw(const string& out_path, const string& in_path, boost::unordered_set<string>& o, const boost::unordered_set<string>& ram_id);
	void create_umm_tmp_from_cluster(const string& out_path, const string& in_path, boost::unordered_set<string>& o, const boost::unordered_map<string, int32_t> cluster_entries, const boost::unordered_set<string>& ram_read_ids);
	void write_non_dup_umm(const string& out_path, const string& in_path, const boost::unordered_map<string, int32_t>& dup_cnt_map);

	// transduction contig
	void run_transduction_contig();
	void create_contig_two_ram(const string& out_path, const string& in_path);
	void create_fa_from_tea_contig(const string& out_path_1, const string& out_path_2, const string& in_path);
	int64_t rfind_the_last_of(int64_t& n_len, const string& str, const char chr);
	int64_t find_the_last_of(int64_t& n_len, const string& str, const char chr);

	void collect_two_ram_map(boost::unordered_map<int64_t, string>& two_ram_map, const string& in_path);
	void collect_two_ram_seq_id_ref_set(boost::unordered_map<string, string>& aligned_map, set<string>& two_ram_seq_id_set, const string& in_path, const boost::unordered_map<int64_t, string>& two_ram_map);
	void collect_aln_sam_repeat(set<string>& repeat_selected_seq_id, boost::unordered_set<string>& repeat_two_ram_id, const string& in_path);
	void collect_aln_ram_sam_repeat(boost::unordered_map<string, string>& repeat_ram_aligned_map, const string& in_path);
	void create_tea_transduction(const string& out_path, set<string>& gold, const boost::unordered_map<string, string>& repeat_ram_aligned_map, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path);

	// orphan
	void run_orphan();
	void collect_gene_sets(set<string>& oo, set<string>& gene, const string& in_path);
	void collect_rescue_sets(set<string>& rescue, set<string>& gene, const string& in_path);
	void create_bed_filtered_intersected(const string& out_path, const set<string>& ooo, const string& in_path);
	void collect_cluster_id_pos_map(boost::unordered_map<string, string>& cluster_id_pos_map, const string& in_path);
	void create_intersected_filtered_insertion(const string& out_path, const boost::unordered_map<string, string> cluster_id_pos_map, const string& in_path);
	void collect_cluster_id_pair_map(boost::unordered_map<string, string>& cluster_id_pair_map, const string& in_path);
	void create_orphan_umm_from_cluster(const string& out_path, const boost::unordered_map<string, string>& cluster_id_pair_map, const string& in_path);

	// orphan contig
	void run_orphan_contig();
	void create_orphan_fa_from_tea_contig(const string& out_path, const string& in_path);
	void create_tea_orphan(const string& out_path, set<string>& gold, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path);

	void clean();

	// post processing
	void post_process();
	void create_orphan_list(const string& out_path, set<string>& orphan, const string& in_path);
	void create_transduction_filtered(const string& out_path, set<string>& trans, const string& in_path, const set<string>& orphan);
	void create_contig_filtered_fa(const string& out_path, const string& in_path);
	void collect_aln_repeat_selected_seq_id(set<string>& repeat_selected_seq_id, const string& in_path);
	void create_tea_contig_tmp(const string& out_path, const string& in_path, const set<string>& o);
	void create_tea_contig_tmp_tmp(const string& out_path, const string& in_path, const set<string>& allset);
	void refine_tea_contig_tmp_tmp(const string& out_path, const string& in_path);
	void create_tea_contig_tmp_tmp_filtered_fa(const string& out_path, const string& in_path);
	void collect_refined_aln_sam_ref(set<string>& rname, boost::unordered_map<string, string>& aligned_map, set<string>& two_ram_seq_id_set, const string& in_path, const boost::unordered_map<int64_t, string>& two_ram_map);
	void collect_refined_aln_sam_repeat(set<string>& rrname, set<string>& oo, set<string>& ooo, const string& in_path, const set<string>& o);
	void create_short_transduction_list(const string& out_path, const set<string>& candidate, const set<string>& o, const set<string>& overlap, const boost::unordered_map<string, string>& ref_aligned_map, const string& in_path);
	void collect_transduction_set(set<string>& o, const string& in_path);
	void create_post_contig_list(const string& out_path, const set<string>& o, const string& in_path);

	void output_raw_file(const string& chr, const string& cl_prefix, const RAMIntervalVector& p_cl, const RAMIntervalVector& n_cl, const multimap<int64_t, int64_t>& pm_cl, const boost::unordered_set<int64_t>& positive_only, const boost::unordered_set<int64_t>& negative_only, const int64_t read_length, const int64_t fragment_size, const bool headless);
	void BAM_to_FASTQ_serial(const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void BAM_to_FASTQ(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void BAM_to_FASTQ__MEM(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void BAM_to_FASTQ__MEM_alt(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);

	void MEMBAM_to_FASTQ(const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void collect_boundaries_un(vector<meerkat::BlockBoundary>& fixed_size_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path, int64_t size_block);
	void collect_boundaries_pos(vector<meerkat::BlockBoundary>& fixed_size_blocks, vector<meerkat::BlockBoundary>& unmapped_included_blocks, vector<meerkat::BlockBoundary>& independent_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path, int64_t size_block);
	void collect_boundaries_alt(vector<BlockOffset>& offset_blocks, const string& a_path, const string& a_bai_path, const string& a_bni_path);
	void create_bni_even_index(const string& a_path, const string& a_bai_path, const string& a_bni_path);
//	void collect_boundaries(vector<meerkat::BlockBoundary>& fixed_size_blocks, const string& input_BAM_name, int64_t size_block, bool verbose = false);


private:
	void _BAM_to_FASTQ(vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void _BAM_to_FASTQ__MEM(vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void _BAM_to_FASTQ__MEM_alt(vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void _BAM_to_FASTQ__MEM_new_logic(vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name, const string& orphan_FASTQ_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void _MEMBAM_to_FASTQ(vector<int64_t>& block_boundary, const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void _MEMBAM_to_FASTQ_serial(const string& input_BAM_name, const string& disc_1_FASTQ_name, const string& disc_2_FASTQ_name);
	void remove_entry_enclosed_with_large_H(vector<BamAlignment>& alns);
	void convert_H_to_S(vector<BamAlignment>& alns);
	void fill_H_to_S(BamAlignment& aln, const AlnSeqQualEntry& aln_seq_entry);
	void split_query_to_segments(vector<string>& seqs, vector<string>& quals, vector<BamAlignment>& alns);
	void format_isize();
	void create_disc_FASTQs();
	void create_um_FASTQs();
	void generate_va_bams();
	void generate_vam_files();
	void generate_um_bams(const string& a_path, const string& a_bai_path, const string& a_bni_path);
	void generate_um_raw_bam_serial();
	void generate_um_raw_bam(vector<meerkat::BlockBoundary>& actual_blocks);
	void remove_duplicates_serial(const string& in_file_name, const string& out_file_name);
	void remove_duplicates(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& out_bam_file_name, const string& out_sam_file_name);
	void generate_umm(const string& in_file_name, const string& out_file_name);
	void write_duplicates(const string& a_path, const string& a_bai_path, const string& a_bni_path, const string& out_file_name);
	void generate_ra_bams();
	void generate_ram_files();
	void generate_ram_file(const string& refbam, const string& rbamf, const string& ramf, const string& disc_1_ra_bam, const string& disc_2_ra_bam, bool exo = false, bool headless = false);
	void load_repeat_mapping(boost::unordered_map<string, string>& h, bool& exo, vector<string>& rabam_files);
	void _load_repeat_mapping(boost::unordered_map<string, string>& h, bool& exo, const int32_t end, const string& rabam);

	// The function, write_ram_and_bam, is not fully implemented, incorrect and inefficient
	void write_ram_and_bam(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo);
	void write_ram_and_bam_serial(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo, bool headless = false);
	void write_ram_and_bam_mem_serial(const string& refbam, const string& rbamf, const string& ramf, const boost::unordered_map<string, string>& h, const bool exo, bool headless = false);

	void generate_cbam_files();
	void _generate_cbam_files_mem();
	void _generate_cbam_files_mem_alt();
	void _generate_cbam_files_mem_alt2();
	void _generate_cbam_files_mem_org();
	void _generate_cbam_files_sampe();
	void generate_cbam_files_serial();
	void is_containing_S_or_H(bool& has_S, bool& has_H, vector<BamAlignment>& alns);
	bool has_S_or_H(const BamAlignment& aln);
	bool has_tag(const BamAlignment& aln, const char tag);

	void find_quality_standard(BamTools::BamReader& a_reader);
	void load_repeat_annotation(boost::unordered_map<string, pair<string, string>>& rannot);
	void load_virus_annotation(map<int64_t, string>& vannot);
	void load_ref_annotation(set<string>& chrl,
			boost::unordered_map<string, pair<string, string>>& rannot,
//			boost::unordered_map<string, boost::unordered_map<string, vector< pair<int64_t, int64_t> > > >& ril_annot,
			boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt,
			boost::unordered_map<string, vector<pair<int64_t, int64_t>>>& gap_annot,
			boost::unordered_map<string, GeneIntervalVector>& gene_annot, bool out_chrl, bool out_gap);
	void load_chr(set<string>& chrl) const;
	void read_gap_rfile(boost::unordered_map<string, vector<pair<int64_t, int64_t>>>& gap_annot, const string& file_name, const set<string>& chrl);
	void read_gene_rfile(boost::unordered_map<string, GeneIntervalVector>& gene_annot, const string& file_name, const set<string>& chrl);
	void load_read_length(boost::unordered_map<string, int32_t>& rl);
	void load_insert_size(boost::unordered_map<string, boost::unordered_map<string, double>>& is, boost::unordered_map<string, int32_t>& rl);
	void load_ram(
			boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram,
			boost::unordered_map<string, pair<string, string>>& rannot,
			const bool rm_dup = false );
	void get_cluster_alt(const string& chr, RAMIntervalVector& cl, vector<RAMRepeatEntry>& sram, boost::unordered_map<string, pair<string, string>>& rannot, const int32_t strand, int64_t gap_cutoff);
	void pair_cluster_alt(multimap<int64_t, int64_t>& pm_cl, RAMIntervalVector& p_cl, RAMIntervalVector& n_cl, const int64_t gap_cutoff, const int64_t read_length, bool stringent_pair = false);

	void count_clipped(
			boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt,
			boost::unordered_map<string, GeneIntervalVector>& gene_annot,
			const string& chr,
			const string& cl_prefix,
			const string& contig_dir,
			const multimap<int64_t, int64_t>& pm_cl,
			const RAMIntervalVector& p_cl,
			const RAMIntervalVector& n_cl,
			boost::unordered_set<int64_t>& positive_only,
			boost::unordered_set<int64_t>& negative_only,
			const int64_t read_length,
			const int64_t fragment_size,
			const int64_t rmasker_filter_margin,
			const int64_t gene_margin,
			const bool headless);

	void count_clipped_append(
			const string& chr,
			const string& cl_prefix,
			const string& contig_dir,
			const int64_t read_length,
			const int64_t fragment_size,
			const bool headless);

	void count_clipped_v(
			boost::unordered_map<string, RefRepeatIntervalVector>& ril_annot_alt,
			map<int64_t, string>& vannot,
			const string& chr,
			const string& cl_prefix,
			const string& contig_dir,
			const multimap<int64_t, int64_t>& pm_cl,
			const RAMIntervalVector& p_cl,
			const RAMIntervalVector& n_cl,
			boost::unordered_set<int64_t>& positive_only,
			boost::unordered_set<int64_t>& negative_only,
			const int64_t read_length,
			const int64_t fragment_size,
			const int64_t rmasker_filter_margin,
			const int64_t gene_margin);

	void get_clipped_entries(vector<ClippedEntry>& clipped_entries, int64_t& max_pos_positive, int64_t& max_pos_negative, int64_t& n_positive_clipped_reads, int64_t& n_negative_clipped_reads, int64_t& n_aligned_clipped_positive, int64_t& n_aligned_clipped_negative, BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const int64_t read_length, const int64_t mid_point);

	void get_clipped_entries_append(vector<ClippedEntry>& clipped_entries, int64_t& max_pos_positive, int64_t& max_pos_negative, int64_t& n_positive_clipped_reads, int64_t& n_negative_clipped_reads, int64_t& n_aligned_clipped_positive, int64_t& n_aligned_clipped_negative, BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const int64_t read_length, const int64_t mid_point);

	void output_clipped_stat(ofstream& out_p_clipped_filename, ofstream& out_n_clipped_filename, ofstream& out_p_mate_rname, ofstream& out_n_mate_rname, ofstream& out_cl, ofstream& out_tea, ofstream& out_clipped, const string& contig_dir, RefRepeatIntervalTree& ref_repeat_interval_tree, RefRepeatIntervalVector& stat_results, GeneIntervalTree& gene_interval_tree, GeneIntervalVector& gene_results, BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const string& prefixed_chr, const int64_t read_length, const int64_t rmasker_filter_margin, const int64_t gene_margin, const int64_t mid_point);

	void output_clipped_stat_append(
			ofstream& out_p_clipped_filename,
			ofstream& out_n_clipped_filename,
			ofstream& out_p_mate_rname,
			ofstream& out_n_mate_rname,
			const string& contig_dir,
			BamTools::BamReader& local_reader,
			const int64_t the_ram_boundary_start,
			const int64_t the_ram_boundary_end,
			const RAMIntervalEntry& positive_entry,
			const RAMIntervalEntry& negative_entry,
			const string& chr,
			const string& prefixed_chr,
			const int64_t read_length,
			const int64_t mid_point);


	void output_clipped_stat_v(ofstream& out_p_clipped_filename, ofstream& out_n_clipped_filename, ofstream& out_p_mate_rname, ofstream& out_n_mate_rname, ofstream& out_cl, ofstream& out_tea, ofstream& out_clipped, const string& contig_dir, RefRepeatIntervalTree& ref_repeat_interval_tree, RefRepeatIntervalVector& stat_results, const map<int64_t, string>& vannot,
				BamTools::BamReader& local_reader, const int64_t the_ram_boundary_start, const int64_t the_ram_boundary_end, const RAMIntervalEntry& positive_entry, const RAMIntervalEntry& negative_entry, const string& chr, const string& prefixed_chr, const int64_t read_length, const int64_t rmasker_filter_margin, const int64_t gene_margin);
	void output_mate_fa(boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram);

	void _output_mate_fa(
			boost::unordered_map<string, vector<string>>& positive_mate_reads,
			boost::unordered_map<string, vector<string>>& negative_mate_reads,
			vector<meerkat::BlockBoundary>& actual_blocks,
			const string& input_BAM_name,
			const multimap<string, AlnPairEntry>& a_positive_repeat_map,
			const multimap<string, AlnPairEntry>& a_negative_repeat_map);

	void _output_mate_fa_ram(
			boost::unordered_map<string, vector<string>>& positive_mate_reads,
			boost::unordered_map<string, vector<string>>& negative_mate_reads,
			const multimap<string, AlnPairEntry>& a_positive_repeat_map,
			const multimap<string, AlnPairEntry>& a_negative_repeat_map,
			boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram);

	void output_mate_fa_v(boost::unordered_map<string, boost::unordered_map<int8_t, vector<RAMRepeatEntry>>>& ram);

	void _output_mate_fa_v(
			boost::unordered_map<string, vector<string>>& positive_mate_reads,
			boost::unordered_map<string, vector<string>>& negative_mate_reads,
			vector<meerkat::BlockBoundary>& actual_blocks, const string& input_BAM_name,
			const boost::unordered_map<string, AlnPairEntry>& a_positive_repeat_map,
			const boost::unordered_map<string, AlnPairEntry>& a_negative_repeat_map);

	void _output_mate_fa_serial(
			boost::unordered_map<string, vector<string>>& positive_mate_reads,
			boost::unordered_map<string, vector<string>>& negative_mate_reads,
			const string& input_BAM_name,
			const boost::unordered_map<string, AlnPairEntry>& a_positive_repeat_map,
			const boost::unordered_map<string, AlnPairEntry>& a_negative_repeat_map);

private:
	int64_t get_number_of_good_qualities(const string& a_qual, const int64_t qcutoff);
	int64_t get_number_of_low_qualities_at_begin(const string& a_qual, const int64_t qcutoff);
	int64_t get_number_of_low_qualities_at_end(const string& a_qual, const int64_t qcutoff);
	int64_t get_cpos(int32_t pos, std::vector<BamTools::CigarOp>& cigar, const string& qual, int8_t strand);
	void get_longest_fa(string& a_contig, const string& fa_path);
	void get_two_longest_fa(pair<string, string>& contigs, const string& fa_path, char type);
	void get_bai_index_path(const string& parent_path, string& bai_path);
	void get_bni_index_path(const string& parent_path, string& bni_path);
private:
	int32_t n_cores;
	int32_t qenc;
	uint8_t min_qual;
	TEAOptionParser options;
	vector<meerkat::BlockBoundary> fixed_size_blocks;
	vector<meerkat::BlockBoundary> independent_blocks;
	vector<meerkat::BlockBoundary> unmapped_included_blocks;
};

// maximal length of a good quality value
inline int64_t TEA::get_number_of_good_qualities(const string& a_qual, const int64_t qcutoff) {
	int64_t n_quals = 0;
	uint32_t limit_val = min_qual + qcutoff;
	int64_t max_quals = 0;
	for(uint64_t qual_id = 0; qual_id < a_qual.size(); ++qual_id) {
		uint8_t a_qual_value = a_qual[qual_id];
		if(a_qual_value > limit_val) {
			++n_quals;
			max_quals = max(n_quals, max_quals);
		} else {
			n_quals = 0;
		}
	}
	return max_quals;
}

inline int64_t TEA::get_number_of_low_qualities_at_begin(const string& a_qual, const int64_t qcutoff) {
	int64_t n_quals = 0;
	uint32_t limit_val = min_qual + qcutoff;
	for(uint64_t qual_id = 0; qual_id < a_qual.size(); ++qual_id) {
		uint8_t a_qual_value = a_qual[qual_id];
		if(a_qual_value == limit_val) {
			++n_quals;
		} else {
			break;
		}
	}
	return n_quals;
}
inline int64_t TEA::get_number_of_low_qualities_at_end(const string& a_qual, const int64_t qcutoff) {
	int64_t n_quals = 0;
	uint32_t limit_val = min_qual + qcutoff;
	int64_t min_id = a_qual.size();
	--min_id;
	for(int64_t qual_id = min_id; qual_id >=0; --qual_id) {
		uint8_t a_qual_value = a_qual[qual_id];
		if(a_qual_value == limit_val) {
			++n_quals;
		} else {
			break;
		}
	}
	return n_quals;
}

inline int64_t TEA::get_cpos(int32_t pos, std::vector<BamTools::CigarOp>& cigar, const string& qual, int8_t strand) {
	int32_t cpos = pos;

	if (1 == strand) {
		return cpos;
	}
	if(cigar.empty()) {
		return cpos;
	}

//	# for the negative strand clipped reads, count the number of base pairs (M & I not D)
	int32_t gap = 0;
	uint64_t c_id = 0;
	if('S' == cigar.front().Type) {
		c_id = 1;
	}
	for(; c_id < cigar.size(); ++c_id) {
		auto& a_cigar = cigar[c_id];
		if('I' == a_cigar.Type || 'S' == a_cigar.Type) {
			continue;
		}
		gap += a_cigar.Length;
	}

	cpos += gap;

	return -cpos;
}

inline void TEA::get_longest_fa(string& a_contig, const string& fa_path) {
	stringstream tmp_lines;
	string line;
	ifstream in(fa_path, ios::binary);
	while(getline(in, line, '\n')) {
		if('>' == line[0]) {
			string the_line = tmp_lines.str();
			if(the_line.size() > a_contig.size()) {
				a_contig = the_line;
			}
			tmp_lines.str(string());
			continue;
		}
		tmp_lines << line;
	}
	string the_line = tmp_lines.str();
	if(the_line.size() > a_contig.size()) {
		a_contig = the_line;
	}
}

inline void TEA::get_two_longest_fa(pair<string, string>& contigs, const string& fa_path, char type) {
	stringstream tmp_lines;
	string line;
	ifstream in(fa_path, ios::binary);
	while(getline(in, line, '\n')) {
		if('>' == line[0]) {
			string the_line = tmp_lines.str() + ":" + type;
			if (the_line.size() > 2){
				if(the_line.size() > contigs.second.size()) {
					contigs.second = the_line;
				}
				if(contigs.second.size() > contigs.first.size()) {
					the_line = contigs.first;
					contigs.first = contigs.second;
					contigs.second = the_line;
				}
			}
			tmp_lines.str(string());
			continue;
		}
		tmp_lines << line;
	}
	string the_line = tmp_lines.str() + ":" + type;
	if (the_line.size() > 2){
		if(the_line.size() > contigs.second.size()) {
			contigs.second = the_line;
		}
		if(contigs.second.size() > contigs.first.size()) {
			the_line = contigs.first;
			contigs.first = contigs.second;
			contigs.second = the_line;
		}
	}
}

inline void TEA::get_bai_index_path(const string& parent_path, string& bai_path) {
	bai_path = parent_path;
	bai_path += ".bai";
	if(boost::filesystem::exists(bai_path)) {
		return;
	}
	string target_path(parent_path);
	target_path.back() = 'i';
	if(boost::filesystem::exists(target_path)) {
		bai_path = target_path;
		return;
	}
	bai_path = options.prefix + ".tmp.bai";
	if (!options.working_dir.empty()) {
		bai_path = options.working_prefix + ".tmp.bai";
	}
	if(boost::filesystem::exists(bai_path)) {
		return;
	}
	string sambamba_cmd = (boost::format("sambamba index -t %d %s %s") % n_cores % parent_path % bai_path).str();
	system(sambamba_cmd.c_str());
}

inline void TEA::get_bni_index_path(const string& parent_path, string& bni_path) {
	string target_path(parent_path);
	target_path[target_path.size() - 2] = 'n';
	target_path[target_path.size() - 1] = 'i';
	if(boost::filesystem::exists(target_path)) {
		bni_path = target_path;
		return;
	}
	bni_path = parent_path;
	bni_path += ".bni";
	if(boost::filesystem::exists(bni_path)) {
		return;
	}

//	bni_path = options.prefix + ".bni";
//	if (!options.working_dir.empty()) {
//		bni_path = options.working_prefix + ".bni";
//	}
}

} /* namespace tea */

#endif /* TEA_TEA_HPP_ */
