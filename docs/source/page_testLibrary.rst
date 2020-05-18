########################################
The tests of HippoUnit in the HBP Validation framework
########################################

{% for test in test_json %}

{{ test["name"] }}
{% for i in range(test["name"]|length) %}-{% endfor %}

.. raw:: html

    <div>
        <style type="text/css">
            .firstColumnBold td:first-child { font-weight: bold; }
        </style>
        <table style="table-layout:fixed; width:100%; word-break: break-all;" class="firstColumnBold">
            <colgroup>
                <col style="width: 25%;" />
                <col>
            </colgroup>
            <tr>
                <td>id</td>
                <td><a href="https://collab.humanbrainproject.eu/#/collab/8123/nav/409289?state=test.{{test['id']}}" target="_blank">{{ test["id"] }}</a></td>
            </tr>
            <tr>
                <td>uri</td>
                <td>{{ test["uri"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>name</td>
                <td>{{ test["name"] }}</td>
            </tr>
            <tr>
                <td>alias</td>
                <td>{{ test["alias"] }}</td>
            </tr>
            <tr>
                <td>author</td>
                <td>
                    {% for name in test["author"] %} {{ name["given_name"] }} {{ name["family_name"] }} {{ ", " if not loop.last }} {% endfor %}
                </td>
            </tr>
            <tr>
                <td>creation_date</td>
                <td>{{ test["creation_date"] }}</td>
            </tr>
            <tr>
                <td>status</td>
                <td><span> {{ test["status"] }}</span></td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>species</td>
                <td>{{ test["species"] }}</td>
            </tr>
            <tr>
                <td>brain_region</td>
                <td>{{ test["brain_region"] }}</td>
            </tr>
            <tr>
                <td>cell_type</td>
                <td>{{ test["cell_type"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>data_location</td>
                <td><a href='{{ test["data_location"] }}' target="_blank">{{ test["data_location"] }}</a></td>
            </tr>
            <tr>
                <td>data_type</td>
                <td>{{ test["data_type"] }}</td>
            </tr>
            <tr>
                <td>data_modality</td>
                <td>{{ test["data_modality"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>test_type</td>
                <td>{{ test["test_type"] }}</td>
            </tr>
            <tr>
                <td>score_type</td>
                <td>{{ test["score_type"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>protocol</td>
                <td>{{ test["protocol"].encode('unicode_escape')|e }}</td>
            </tr>
        </table>
    </div>

    {% for test_instance in test["codes"] %}

    <div>
        <style type="text/css">
            .firstColumnBold td:first-child { font-weight: bold; }
        </style>
        <table style="table-layout:fixed; width:100%; word-break: break-all;" class="firstColumnBold">
            <colgroup>
                <col style="width: 25%;" />
                <col>
            </colgroup>
            <tr class="card-panel orange lighten-4">
            <th style="text-align:center" colspan="2">Test Instance: {{ test_instance["version"] }}</th>
            </tr>
            <tr>
                <td>id</td>
                <td><a href="https://collab.humanbrainproject.eu/#/collab/8123/nav/409289?state=test.{{test['id']}}">{{ test_instance["id"] }}</a></td>
            </tr>
            <tr>
                <td>uri</td>
                <td>{{ test_instance["uri"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>version</td>
                <td>{{ test_instance["version"] }}</td>
            </tr>
            <tr>
                <td>repository</td>
                <td><a href="{{ test_instance["repository"] }}">{{ test_instance["source"] }}</a></td>
            </tr>
            <tr>
                <td>path</td>
                <td><a href="./api/{{ test_instance["path"] }}.html#{{ test_instance["path"] }}">{{ test_instance["path"] }}</a></td>
            </tr>
            <tr>
                <td>timestamp</td>
                <td>{{ test_instance["timestamp"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>parameters</td>
                <td>{{ test_instance["parameters"] }}</td>
            </tr>
            <tr>
                <td>description</td>
                <td>{{ test_instance["description"] }}</td>
            </tr>
        </table>
    </div>

    {% endfor %}

{% endfor %}