#!/usr/bin/env bash

master_list="http://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=fungi"

###############################################################################

function login() {
  local uname=""; # ADD JGI USERNAME/EMAIL
  local pass="";  # ADD JGI PASSWORD
  curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode "login=$uname" --data-urlencode "password=$pass" -c cookies > /dev/null
}

###############################################################################

function get_published_list() {
  ./select_published_jgi_genomes.sh
}

###############################################################################

function filter_published() {
  local xml="$1"
  local list="$2"
  local out="$3"

  grep -i -f "$list" "$xml" > "$out"
}


###############################################################################

function get_master_list() {
  local xml="$1"
  curl 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism=fungi' -b cookies \
   | xmllint --xpath '//folder[@name="Assembled scaffolds (unmasked)"] | //folder[@name="Genes"]'
   > $xml
}

###############################################################################

function get_relevant_entries() {

  local xml="$1"
  local out="$2"

   xmllint --xpath '//folder[@name="Assembled scaffolds (unmasked)"] | //folder[@name="Genes"]' "$xml" > "$out"

}

###############################################################################

function get_shortname() {
  local data="$1"

  echo "<l>$data</l>" | xmllint --xpath "//l/file/@url" - | sed -e 's/" /\n/g' | rev | cut -d= -f1 | rev | tr -d '"' | cut -d/ -f2
}

###############################################################################

function get_longname() {
  local data="$1"

  echo "<l>$data</l>" | xmllint --xpath "//l/file/@label" - | sed -e 's/" /\n/g' | cut -d= -f1 --complement | tr -d '"'
}

###############################################################################

function get_timestamp() {
  local data="$1"

  date --date="`echo \"<l>$data</l>\" | xmllint --xpath '//l/file/@timestamp' - | tr -d '"' | cut -d= -f2`" +%s
}

###############################################################################

function get_url() {
  local data="$1"

  echo "<l>$data</l>" | xmllint --xpath "//l/file/@url" - | cut -d= -f1 --complement | tr -d '"'
}

###############################################################################

function check() {
  local url="$1"

  if [ ! -z `echo $url | grep -i -e "[._]ests\?[._]\?"` ]; then
    echo "1"
    return 0
  fi

  extensions=(`echo $url | rev | cut -d/ -f1 | rev | cut -d. -f1 --complement | tr '.' '\n' | tac | tr '\n' ' '`)

  for ext in ${extensions[@]}; do
    case $ext in
      gff|fa|fasta) echo "0"; return 0;;
      gff3) echo "0"; return 0;;
      *)             ;;
    esac
  done
  echo "not_ok"
  return 0
}

###############################################################################

function number_target() {
  local target="$1"
  local number="$2"

  if [ -z "$number" ]; then
  
    if [ ! -e "$target" ]; then
      echo "$target"
    else
      number_target "$target" 1
    fi
  else
    if [ ! -e "$target.$number" ]; then
      echo "$target.$number"
    else
      number_target "$target" $((number+1))
    fi
  fi
  
}

###############################################################################

function download_file() {

  local file=`echo "$1" | sed -e 's/&amp;/\&/g'`
  local target="$2"

  mkdir -p `dirname "$target"`
  curl -s "http://genome.jgi.doe.gov/$file" -b cookies > $target

  return $?
}

###############################################################################
###############################################################################
###############################################################################

xml="master.xml"
fxml="master_filt.xml"

#login
#get_master_list "$xml"
get_relevant_entries "$xml" "$xml.rel"
#get_published_list

filter_published "$xml.rel" published_genomes_names "$fxml"

echo -e "#long_name\tshort_name\ttimestamp\tfile" > $(pwd)/download_list
while read line; do
  short_name=`get_shortname "$line"`
  long_name=`get_longname "$line"`
  timestamp=`get_timestamp "$line"`
  url=`get_url "$line"`
 
  if [ "`check "$url"`" == "0" ]; then
    target="$(pwd)/data/`echo $long_name | tr ' ' '_'`/$short_name/`basename $url`"
    target_num=`number_target "$target"`
    download_file $url $target_num
    if [ $? -eq 0 ]; then
      echo -e "\"$long_name\"\t$short_name\t$timestamp\t$target_num" >> $(pwd)/download_list
      echo -e "\"$long_name\"\t$short_name\t$timestamp\t$target_num"
    else
      echo -e "\"$long_name\"\t$short_name\t$timestamp\t$target_num" >> $(pwd)/download_list.failed
    fi
  fi
done < <(cat $fxml | grep -v all | grep -v modifications | grep -v '^#')
