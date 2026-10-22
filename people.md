---
title: People
head_title: people
layout: page
permalink: /people/
---


<div class="people-grid">
{% for person in site.data.people %}
	{% include person.html person=person %}
{% endfor %}
</div>

