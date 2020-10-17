package uniprot

import (
	"fmt"
	"net/url"
	"strings"

	"github.com/tikz/bio/http"
)

// SearchQuery fetches the search results for a given query and returns a list of UniProt IDs.
func SearchQuery(query string) (ids []string, err error) {
	url := "https://www.uniprot.org/uniprot/?query=" + url.QueryEscape(query) + "&format=list"
	fmt.Println(url)
	raw, err := http.Get(url)
	if err != nil {
		return nil, fmt.Errorf("UniProt search query %v: %v", query, err)
	}

	return strings.Split(strings.TrimSpace(string(raw)), "\n"), err
}
