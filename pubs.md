---
title: Publications
head_title: pubs
layout: page
permalink: /pubs/
---

<ol class="publication-list" reversed>
{% for publication in site.data.publications %}
  <li class="publication">
    <strong class="publication-title">{{ publication.title }}</strong>
    <div class="publication-details">
      <span class="publication-venue"><em>{{ publication.journal }}</em>, {{ publication.year }}</span>
      {% for link in publication.links %}
        <a href="{{ link.url }}">[{{ link.label }}]</a>
      {% endfor %}
    </div>
  </li>
{% endfor %}
</ol>
